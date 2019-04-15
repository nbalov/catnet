/*
 *  catnet : categorical Bayesian network inference
 *  Copyright (C) 2009--2010  Nikolay Balov
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  http://www.gnu.org/licenses/gpl-2.0.html
 */

/*
 * rcatnet.cpp
 *
 *  Created on: Sep 21, 2009
 *      Author: nbalov
 */

/* 
 * version 1.15.1  12dec2016
 */

#include "utils.h"
#include "rcatnet_sa.h"
#include "rcatnet.h"

using namespace std;

THREAD_PROC_DEFINE(CatnetSearchSaThreadProc, pParam) {
	SEARCH_PARAMETERS *pSearchParams = (SEARCH_PARAMETERS*)pParam;
	if(pSearchParams && pSearchParams->m_pCaller) {
		CATNET_SEARCH2<char, MAX_NODE_NAME, double> *pCaller = 
			(CATNET_SEARCH2<char, MAX_NODE_NAME, double>*)pSearchParams->m_pCaller;
		pCaller->estimate(pSearchParams);
		pCaller->_exit_thread((void*)0);
	}
	return(0);
}

double *catnetPairwiseCondLikelihood(SEXP rSamples, SEXP rPerturbations);

RCatnetSearchSA::RCatnetSearchSA() {
	m_nDrives = 0;
	m_pDrives = NULL;
	m_pTestOrder = NULL;
	m_pTestOrderInverse = NULL;
	m_bUseCache = 1;
	m_pSearchParams = NULL;
	m_pOptOrder = NULL;
	m_nOptNets = 0;
	m_pOptNets = NULL;
}

RCatnetSearchSA::~RCatnetSearchSA() {
	_release();
}

void RCatnetSearchSA::_release() {
	int i;
	
	if (m_pTestOrder) {
		for (i = 0; i < m_nDrives; i++) {
			if (m_pTestOrder[i])
				CATNET_FREE( m_pTestOrder[i]);
		}
		CATNET_FREE( m_pTestOrder);
		m_pTestOrder = NULL;
	}

	if (m_pTestOrderInverse) {
		for (i = 0; i < m_nDrives; i++) {
			if (m_pTestOrderInverse[i])
				CATNET_FREE( m_pTestOrderInverse[i]);
		}
		CATNET_FREE( m_pTestOrderInverse);
		m_pTestOrderInverse = NULL;
	}

	if (m_pSearchParams) {
		for (i = 0; i < m_nDrives; i++)
			if (m_pSearchParams[i])
				delete m_pSearchParams[i];
		CATNET_FREE( m_pSearchParams);
		m_pSearchParams = NULL;
	}

	if (m_pDrives) {
		for (i = 0; i < m_nDrives; i++)
			if (m_pDrives[i])
				delete m_pDrives[i];
		CATNET_FREE( m_pDrives);
		m_pDrives = NULL;

	}
	m_nDrives = 0;

	if (m_pOptOrder)
		CATNET_FREE( m_pOptOrder);
	m_pOptOrder = NULL;

	if (m_pOptNets && m_nOptNets > 0) {
		for (i = 0; i < m_nOptNets; i++) {
			if (m_pOptNets[i]) {
				delete m_pOptNets[i];
			}
			m_pOptNets[i] = NULL;
		}
		CATNET_FREE( m_pOptNets);
		m_pOptNets = NULL;
		m_nOptNets = 0;
	}

}

int *RCatnetSearchSA::_genOrder(const int *porder, int norder, int shuffles, int bjump) {
	int i, *neworder, sh, n1, n2;
	double u;
	if (norder < 1)
		return 0;
	neworder = (int*) CATNET_MALLOC(norder * sizeof(int));
	if (!neworder)
		return 0;
	if (shuffles <= 0) {
		_gen_permutation<int> (neworder, norder);
		return neworder;
	}
	if (!porder || shuffles < 0)
		return 0;
	int *torder = (int*) CATNET_MALLOC(norder * sizeof(int));
	if (!torder) {
		CATNET_FREE(neworder);
		return 0;
	}
	memcpy(neworder, porder, norder * sizeof(int));
	GetRNGstate();
	for (sh = 0; sh < shuffles; sh++) {
		memcpy(torder, neworder, norder * sizeof(int));
		u = (double) unif_rand();
		n1 = (int) (u * norder);
		if (bjump) {
			n2 = n1;
			while (n2 == n1) {
				u = (double) unif_rand();
				n2 = (int) (u * norder);
			}
		} else {
			n2 = n1 + 1;
			if (n1 >= norder - 1)
				n2 = 0;
		}
		if (n1 < n2) {
			if (n1 > 0)
				for (i = 0; i < n1; i++)
					neworder[i] = torder[i];
			for (i = n1; i < n2; i++)
				neworder[i] = torder[i + 1];
			neworder[n2] = torder[n1];
			if (n2 < norder - 1)
				for (i = n2 + 1; i < norder; i++)
					neworder[i] = torder[i];
		} else {
			if (n2 > 0)
				for (i = 0; i < n2; i++)
					neworder[i] = torder[i];
			for (i = n2; i < n1; i++)
				neworder[i + 1] = torder[i];
			neworder[n2] = torder[n1];
			if (n1 < norder - 1)
				for (i = n1 + 1; i < norder; i++)
					neworder[i] = torder[i];
		}
	}
	PutRNGstate();
	CATNET_FREE(torder);
	return neworder;
}

int *RCatnetSearchSA::_genOrderFormDirProbs(const int *porder, int numnodes, double *matEdgeLiks, double *pOrderProb) {
	int i, j, k, *neworder, *flags;
	double *probs, faux, fsum;
	if (numnodes < 1 || !matEdgeLiks || !pOrderProb)
		return 0;
	neworder = (int*) CATNET_MALLOC(numnodes*sizeof(int));
	if (!neworder)
		return 0;
	memset(neworder, 0, numnodes*sizeof(int));
	flags = (int*) CATNET_MALLOC(numnodes*sizeof(int));
	if (!flags) {
		CATNET_FREE(neworder);
		return 0;
	}
	probs = (double*) CATNET_MALLOC(numnodes*sizeof(double));
	if (!probs) {
		CATNET_FREE(flags);
		CATNET_FREE(neworder);
		return 0;
	}

	*pOrderProb = 1;
	neworder[0] = 0;
	GetRNGstate();
	for(k = 1; k < numnodes; k++) {
		fsum = 0;
		for(i = 0; i <= k; i++) {
			faux = 1;
			if(i > 0) {
				for(j = 0; j < i; j++) 
					faux *= matEdgeLiks[neworder[j]*numnodes+k];
			}
			if(i < k) {
				for(j = i; j < k; j++) 
					faux *= matEdgeLiks[k*numnodes+neworder[j]];
			}
			probs[i] = faux;
			fsum += faux;
		}
		faux = (double)fsum * (double) unif_rand();
		fsum = 0;
		for(i = 0; i < k; i++) {
			fsum += probs[i];
			if(fsum >= faux)
				break;
		}
		*pOrderProb *= probs[i];
		if(i > 0)
			memcpy(flags, neworder, i*sizeof(int));
		flags[i] = k;
		if(i < k)
			memcpy(flags + i + 1, neworder + i, (k-i)*sizeof(int));
		memcpy(neworder, flags, (k+1)*sizeof(int));
	}
	PutRNGstate();
	// order should be with indices in [1, numnodes]
	for(i = 0; i < numnodes; i++)
		neworder[i]++;
	if(*pOrderProb > 0)
		*pOrderProb = (double)log((double)*pOrderProb);
	else
		*pOrderProb = (double)-FLT_MAX;
	CATNET_FREE(flags);
	CATNET_FREE(probs);
	return neworder;
}

SEXP RCatnetSearchSA::search(SEXP rNodeNames, SEXP rSamples,
		SEXP rPerturbations, SEXP rMaxParents, SEXP rParentSizes,
		SEXP rMaxComplexity, SEXP rNodeCats, SEXP rParentsPool, SEXP rFixedParentsPool, SEXP rMaxParentsPool, 
		SEXP rMatEdgeLiks, SEXP rDirProbs, 
		SEXP rModel, SEXP rStartOrder, SEXP rTempStart, SEXP rTempCoolFact,
		SEXP rTempCheckOrders, SEXP rMaxIter, SEXP rOrderShuffles,
		SEXP rStopDiff, SEXP rThreads, SEXP rUseCache, SEXP rEcho) {

	int n, nn, i, j, k, len, maxComplexity, numnets, inet, echo, maxParentsPool;
	int *pRsamples, *pRperturbations, *pSamples, *pPerturbations,
			**parentsPool, **fixedParentsPool, *pPool, *pParentSizes;
	double complx, fLogLik, fOptLogLik, deltaLogLik, ftemp;
	int tempCheckOrders, minShuffles, bjump, niter, maxIter, nstop,
			nstopCounter, naccept, stepiters, stepaccept, nLastChangedTemp;
	double tempCur, tempStart, tempCoolFact, orderShuffles, stopDiff, 
			acceptprob;

	MUTEX m_cache_mutex;

	int nModelSel; // 0-maxcomplx, 1-AIC, 2-BIC
	int nCatnets, nCurNet;
	CATNET<char, MAX_NODE_NAME, double> **pCatnets, *pCurNet, *pOptNet;
	const char *pstr;
	char **pNodeNames;

	double *matEdgeLiks, *pMatEdgeLiks, ordProb, *newOrdProb, *pDirProbs, *matNodeCondLiks, *pNodeCondLiks;
	
	RCatnet rcatnet;
	SEXP dim, rnodecat, rparpool, cnetlist, cnetnode;

	_release();

	if (!isMatrix(rSamples))
		error("Data is not a matrix");

	PROTECT(rMaxParents = AS_INTEGER(rMaxParents));
	PROTECT(rMaxComplexity = AS_INTEGER(rMaxComplexity));
	PROTECT(rMaxParentsPool = AS_INTEGER(rMaxParentsPool));
	PROTECT(rTempStart = AS_NUMERIC(rTempStart));
	PROTECT(rTempCoolFact = AS_NUMERIC(rTempCoolFact));
	PROTECT(rTempCheckOrders = AS_INTEGER(rTempCheckOrders));
	PROTECT(rMaxIter = AS_INTEGER(rMaxIter));
	PROTECT(rOrderShuffles = AS_NUMERIC(rOrderShuffles));
	PROTECT(rStopDiff = AS_NUMERIC(rStopDiff));
	PROTECT(rThreads = AS_INTEGER(rThreads));
	PROTECT(rUseCache = AS_LOGICAL(rUseCache));
	PROTECT(rModel = AS_CHARACTER(rModel));
	PROTECT(rEcho = AS_LOGICAL(rEcho));

	m_maxParentSet = INTEGER_POINTER(rMaxParents)[0];
	maxComplexity = INTEGER_POINTER(rMaxComplexity)[0];
	maxParentsPool = INTEGER_POINTER(rMaxParentsPool)[0];
	tempStart = NUMERIC_POINTER(rTempStart)[0];
	tempCoolFact = NUMERIC_POINTER(rTempCoolFact)[0];
	tempCheckOrders = INTEGER_POINTER(rTempCheckOrders)[0];
	maxIter = INTEGER_POINTER(rMaxIter)[0];
	orderShuffles = NUMERIC_POINTER(rOrderShuffles)[0];
	stopDiff = NUMERIC_POINTER(rStopDiff)[0];
	m_nDrives = INTEGER_POINTER(rThreads)[0];
	if (m_nDrives < 1)
		m_nDrives = 1;
	m_bUseCache = LOGICAL(rUseCache)[0];
	if(m_bUseCache && m_nDrives > 8) {
		// with many threads in parallel the cache is presumably inefficient - disable it!
		m_bUseCache = 0;
		warning("Cache is disabled, too many threads");
	}
	
	echo = LOGICAL(rEcho)[0];

	nModelSel = 0;
	if (!strncmp(CHARACTER_VALUE(rModel), "AIC", 3))
		nModelSel = 1;
	if (!strncmp(CHARACTER_VALUE(rModel), "BIC", 3))
		nModelSel = 2;

	UNPROTECT(13);

	PROTECT(rSamples = AS_INTEGER(rSamples));
	pRsamples = INTEGER(rSamples);

	dim = GET_DIM(rSamples);
	m_numNodes = INTEGER(dim)[0];
	m_numSamples = INTEGER(dim)[1];

	m_pOptOrder = (int*) CATNET_MALLOC(m_numNodes * sizeof(int));
	PROTECT(rStartOrder = AS_INTEGER(rStartOrder));
	if (length(rStartOrder) < m_numNodes) {
		warning("Invalid nodeStartOrder parameter - reset to default node order.");
		for (i = 0; i < m_numNodes; i++)
			m_pOptOrder[i] = i + 1;
	} else {
		if (m_pOptOrder) {
			memcpy(m_pOptOrder, INTEGER(rStartOrder), m_numNodes * sizeof(int));
			for (i = 0; i < m_numNodes; i++) {
				if(m_pOptOrder[i] <= 0 || m_pOptOrder[i] > m_numNodes) {
					error("Invalid startOrder parameter");
				}
			}
		}
	}
	UNPROTECT(1); /* rStartOrder */

	pNodeNames = NULL;
	if (!isNull(rNodeNames)) {
		PROTECT(rNodeNames = AS_VECTOR(rNodeNames));
		if (length(rNodeNames) == m_numNodes) {
			pNodeNames = (char**) CATNET_MALLOC(m_numNodes * sizeof(char*));
			for (i = 0; i < m_numNodes; i++) {
				pstr = CHARACTER_VALUE(VECTOR_ELT(rNodeNames, i));
				pNodeNames[i] = (char*) CATNET_MALLOC((strlen(pstr) + 1) * sizeof(char));
				if (pNodeNames[i] && pstr)
					memcpy(pNodeNames[i], pstr, (strlen(pstr) + 1) * sizeof(char));
			}
		}
		UNPROTECT(1); /* rNodeNames */
	}

	bjump = 0;
	if (orderShuffles < 0) {
		bjump = 1;
		orderShuffles = -orderShuffles;
	}
	minShuffles = (int) (orderShuffles);
	orderShuffles = orderShuffles - minShuffles;

	m_pTestOrder        = (int**) CATNET_MALLOC(m_nDrives * sizeof(int*));		
	m_pTestOrderInverse = (int**) CATNET_MALLOC(m_nDrives * sizeof(int*));
	newOrdProb          = (double*) CATNET_MALLOC(m_nDrives * sizeof(double));

	m_pDrives = (CATNET_SEARCH2<char, MAX_NODE_NAME, double>**) CATNET_MALLOC(
			m_nDrives * sizeof(CATNET_SEARCH2<char, MAX_NODE_NAME, double>*));

	if (!m_pDrives || !m_pTestOrder || !m_pTestOrderInverse || !newOrdProb) {
		if (m_pTestOrder) 
			CATNET_FREE(m_pTestOrder);
		m_pTestOrder = 0;
		if (m_pTestOrderInverse) 
			CATNET_FREE(m_pTestOrderInverse);
		m_pTestOrderInverse = 0;
		if (newOrdProb) 
			CATNET_FREE(newOrdProb);
		if (m_pDrives)
			CATNET_FREE(m_pDrives);
		m_pDrives = 0;
		UNPROTECT(1); // rSamples
		return R_NilValue;
	}
	
	memset(newOrdProb, 0, m_nDrives * sizeof(double));

	memset(m_pTestOrder, 0, m_nDrives * sizeof(int*));

	for (n = 0; n < m_nDrives; n++) {
		m_pTestOrderInverse[n] = (int*) CATNET_MALLOC(m_numNodes * sizeof(int));
	}

	for (n = 0; n < m_nDrives; n++) {
		m_pDrives[n] = new CATNET_SEARCH2<char, MAX_NODE_NAME, double> ();
	}

	if (m_bUseCache) {
		MUTEX_INIT(m_cache_mutex);
	}

	if (!isNull(rParentSizes) && length(rParentSizes) == m_numNodes)
		PROTECT(rParentSizes = AS_INTEGER(rParentSizes));
	if (!isNull(rPerturbations))
		PROTECT(rPerturbations = AS_INTEGER(rPerturbations));
	if(!isNull(rNodeCats))
		PROTECT(rNodeCats = AS_LIST(rNodeCats));
	if (!isNull(rParentsPool) && length(rParentsPool) == m_numNodes)
		PROTECT(rParentsPool = AS_LIST(rParentsPool));
	if (!isNull(rFixedParentsPool) && length(rFixedParentsPool) == m_numNodes)
		PROTECT(rFixedParentsPool = AS_LIST(rFixedParentsPool));
	if (!isNull(rMatEdgeLiks) && length(rMatEdgeLiks) == m_numNodes*m_numNodes)
		PROTECT(rMatEdgeLiks = AS_NUMERIC(rMatEdgeLiks));
	if (!isNull(rDirProbs) && length(rDirProbs) == m_numNodes*m_numNodes)
		PROTECT(rDirProbs = AS_NUMERIC(rDirProbs));

	pNodeCondLiks = 0;
	if(maxParentsPool >= 1) {
		pNodeCondLiks = catnetPairwiseCondLikelihood(rSamples, rPerturbations);
		if(m_bUseCache && 4*maxParentsPool < m_numNodes) {
			m_bUseCache = 0;
			warning("Cache is disabled, maxParentsPool<<numNodes");
		}
	}

	m_pSearchParams = (SEARCH_PARAMETERS**) CATNET_MALLOC(m_nDrives * sizeof(SEARCH_PARAMETERS*));
	for (n = 0; n < m_nDrives; n++) {
		m_pSearchParams[n] = new SEARCH_PARAMETERS(m_numNodes, m_numSamples,
				m_maxParentSet, maxComplexity, (m_nDrives < 2) ? echo : 0, 
				!isNull(rNodeCats), 
				!isNull(rParentSizes), !isNull(rPerturbations), 
				!isNull(rParentsPool), !isNull(rFixedParentsPool), 
				!isNull(rMatEdgeLiks), maxParentsPool, 
				m_bUseCache ? &m_cache_mutex : NULL, m_pDrives[n]);
	}

	fOptLogLik = -FLT_MAX;
	m_nOptNets = 0;
	m_pOptNets = NULL;
	pOptNet = NULL;

	nstop = 0;
	nstopCounter = 0;
	naccept = 0;
	tempCur = tempStart;
	nLastChangedTemp = 0;

	if (echo && m_numNodes < 32) {
		Rprintf("Start order: ");
		for (i = 0; i < m_numNodes; i++)
			Rprintf("%d ", m_pOptOrder[i]);
		Rprintf("\n");
	}

	ordProb = 0;
	niter = 0;
	while (niter < maxIter) {

		for (n = 0; n < m_nDrives; n++) {
			if (m_pTestOrder[n])
				CATNET_FREE( m_pTestOrder[n]);

			if (!isNull(rDirProbs) && length(rDirProbs) == m_numNodes*m_numNodes) {
				pDirProbs = REAL(rDirProbs);
				m_pTestOrder[n] = _genOrderFormDirProbs((const int*)m_pOptOrder, m_numNodes, pDirProbs, newOrdProb + n);
				
			}
			else {
				m_pTestOrder[n] = _genOrder((const int*)m_pOptOrder, m_numNodes, minShuffles
					+ _gen_binomial(1, orderShuffles), bjump);
			}
			if (!m_pTestOrder[n])
				error("_genOrder returns an error");

			nstop = 0;
			if (n > 0) {
				for (nn = 0; nn < n; nn++) {
					nstop = 1;
					for (i = 0; i < m_numNodes; i++) {
						if (m_pTestOrder[nn][i] != m_pTestOrder[n][i]) {
							nstop = 0;
							break;
						}
					}
					if (nstop)
						break;
				}
			}
			if (nstop) 
				continue;

			for (i = 0; i < m_numNodes; i++)
				m_pTestOrderInverse[n][m_pTestOrder[n][i] - 1] = i + 1;

			if (m_bUseCache)
				m_pDrives[n]->setCacheParams(m_numNodes, m_maxParentSet,
						m_pTestOrder[n], m_pTestOrderInverse[n]);

			if (!isNull(rParentSizes) && length(rParentSizes) == m_numNodes) {
				pParentSizes = m_pSearchParams[n]->m_pParentSizes;
				if (pParentSizes && INTEGER(rParentSizes))
					memcpy(pParentSizes, INTEGER(rParentSizes), m_numNodes
						* sizeof(int));
				for (i = 0; i < m_numNodes; i++)
					pParentSizes[i] = INTEGER(rParentSizes)[m_pTestOrder[n][i]-1];
			}

			pSamples = m_pSearchParams[n]->m_pSamples;
			for (j = 0; j < m_numSamples; j++) {
				for (i = 0; i < m_numNodes; i++) {
					pSamples[j * m_numNodes + i] = pRsamples[j * m_numNodes
							+ m_pTestOrder[n][i] - 1];
					if(R_IsNA(pSamples[j*m_numNodes + i]) || pSamples[j*m_numNodes + i] < 1)
						pSamples[j*m_numNodes + i] = CATNET_NAN;
				}
			}

			if (!isNull(rPerturbations)) {
				pPerturbations = m_pSearchParams[n]->m_pPerturbations;
				pRperturbations = INTEGER(rPerturbations);
				for (j = 0; j < m_numSamples; j++) {
					for (i = 0; i < m_numNodes; i++)
						pPerturbations[j * m_numNodes + i] = pRperturbations[j
								* m_numNodes + m_pTestOrder[n][i] - 1];
				}
			}

			if(!isNull(rNodeCats)) {
				for(i = 0; i < m_numNodes; i++) {
					rnodecat = AS_INTEGER(VECTOR_ELT(rNodeCats, (int)(m_pTestOrder[n][i] - 1)));
					len = length(rnodecat);
					if(isVector(rnodecat) && len > 0) {
						m_pSearchParams[n]->m_pNodeNumCats[i] = len;
						m_pSearchParams[n]->m_pNodeCats[i] = (int*)CATNET_MALLOC(len*sizeof(int));
						for(j = 0; j < len; j++)
							m_pSearchParams[n]->m_pNodeCats[i][j] = INTEGER(rnodecat)[j];
					}
				}
			}

			if (!isNull(rParentsPool) && length(rParentsPool) == m_numNodes) {
				parentsPool = m_pSearchParams[n]->m_parentsPool;

				for (i = 0; i < m_numNodes; i++) {
					rparpool = AS_INTEGER(VECTOR_ELT(rParentsPool,
							(int) (m_pTestOrder[n][i] - 1)));
					len = length(rparpool);
					if (isVector(rparpool) && len > 0 && len <= m_numNodes) {
						pPool = INTEGER(rparpool);
						for (j = 0; j < len; j++) {
							if (pPool[j] > 0 && pPool[j] <= m_numNodes) {
								for (k = 0; k < m_numNodes; k++)
									if (pPool[j] == m_pTestOrder[n][k])
										break;
								if (k < m_numNodes)
									parentsPool[i][j] = k;
								else
									parentsPool[i][j] = -1;
							} else
								parentsPool[i][j] = -1;
						}
						for (; j < m_numNodes; j++)
							parentsPool[i][j] = -1;
					} else {
						//warning("Invalid parentsPool parameter");
						for (j = 0; j < m_numNodes; j++)
							parentsPool[i][j] = -1;
						parentsPool[i][0] = i;
					}
				}
			}

			if (!isNull(rFixedParentsPool) && length(rFixedParentsPool) == m_numNodes) {
				fixedParentsPool = m_pSearchParams[n]->m_fixedParentsPool;
				for (i = 0; i < m_numNodes; i++) {
					rparpool = AS_INTEGER(VECTOR_ELT(rFixedParentsPool,
							(int) (m_pTestOrder[n][i] - 1)));
					len = length(rparpool);
					if (isVector(rparpool) && len > 0 && len <= m_numNodes) {
						if (m_maxParentSet < len)
							m_maxParentSet = len;
						pPool = INTEGER(rparpool);
						for (j = 0; j < len; j++) {
							if (pPool[j] > 0 && pPool[j] <= m_numNodes) {
								for (k = 0; k < m_numNodes; k++)
									if (pPool[j] == m_pTestOrder[n][k])
										break;
								if (k < m_numNodes)
									fixedParentsPool[i][j] = k;
								else
									fixedParentsPool[i][j] = -1;
							} else
								fixedParentsPool[i][j] = -1;
						}
						for (; j < m_numNodes; j++)
							fixedParentsPool[i][j] = -1;
					} else {
						for (j = 0; j < m_numNodes; j++)
							fixedParentsPool[i][j] = -1;
					}
				}
			}

			if (!isNull(rMatEdgeLiks)) {
				matEdgeLiks = m_pSearchParams[n]->m_matEdgeLiks;
				pMatEdgeLiks = REAL(rMatEdgeLiks);
				for(j = 0; j < m_numNodes; j++) {
					for(i = 0; i < m_numNodes; i++) {
						matEdgeLiks[j*m_numNodes + i] = pMatEdgeLiks[(m_pTestOrder[n][j] - 1)*m_numNodes + m_pTestOrder[n][i] - 1];
					}
				}
			}

			if (maxParentsPool >= 1 && pNodeCondLiks) {
				matNodeCondLiks = m_pSearchParams[n]->m_matNodeCondLiks;
				for(j = 0; j < m_numNodes; j++) {
					for(i = 0; i < m_numNodes; i++) {
						matNodeCondLiks[j*m_numNodes + i] = pNodeCondLiks[(m_pTestOrder[n][j] - 1)*m_numNodes + m_pTestOrder[n][i] - 1];
					}
				}
			}

			m_pDrives[n] -> _start_thread(CatnetSearchSaThreadProc, m_pSearchParams[n]);
		} // for(n = 0; n < m_nDrives; n++)

		nstop = 1;
		stepiters = m_nDrives;
		stepaccept = 0;

		for (n = 0; n < m_nDrives; n++) {

			if (!m_pDrives[n] ->_is_running())
				continue;
			
			int threadres = m_pDrives[n] -> _join_thread();
			threadres = m_pDrives[n] -> _stop_thread();

			if(stepaccept > 0) {
				// if has acceptance by previous drive, discard the rest
				//printf("break at %d\n", stepiters);
				continue; // must not break it, wait for the others to finish
			}

			nstop = 0;

			nCatnets = 0;
			if (threadres == 0) {
				nCatnets = m_pDrives[n]->numCatnets();
				pCatnets = m_pDrives[n]->catnets();
				//printf("nCatnets = %d,  %p\n", nCatnets, pCatnets);
			}
			if (nCatnets > 0 && pCatnets) {

				pCurNet = NULL;
				switch (nModelSel) {
				case 1: // AIC
					fLogLik = -FLT_MAX;
					for (nCurNet = 0; nCurNet < nCatnets; nCurNet++) {
						if (!pCatnets[nCurNet])
							continue;
						// net complexity = nCurNet + m_numNodes*(1 + maxCategories)
						complx = m_numSamples*pCatnets[nCurNet]->loglik()-pCatnets[nCurNet]->complexity();
						if (complx > fLogLik) {
							fLogLik = complx;
							pCurNet = pCatnets[nCurNet];
						}
					}
					break;
				case 2: // BIC
					fLogLik = -FLT_MAX;
					ftemp = 0.5 * (double)log((double)m_numSamples);
					for (nCurNet = 0; nCurNet < nCatnets; nCurNet++) {
						if (!pCatnets[nCurNet])
							continue;
						complx = m_numSamples*pCatnets[nCurNet]->loglik()-ftemp*pCatnets[nCurNet]->complexity();
						if (complx > fLogLik) {
							fLogLik = complx;
							pCurNet = pCatnets[nCurNet];
						}
					}
					break;
				}

				if(!pCurNet) {
					// optimize w.r.t. max complexity network
					pCurNet = NULL;
					for (nCurNet = nCatnets - 1; nCurNet >= 0; nCurNet--)
						if (pCatnets[nCurNet]) {
							pCurNet = pCatnets[nCurNet];
							break;
						}
				}
				// if nothing has been found
				if(!pCurNet)
					continue;

				fLogLik = pCurNet->loglik();
				if (!pOptNet) {
					pOptNet = pCurNet;
					fOptLogLik = fLogLik;
					if (m_pOptNets) {
						for (i = 0; i < m_nOptNets; i++)
							if (m_pOptNets[i]) {
								delete m_pOptNets[i];
							}
						CATNET_FREE( m_pOptNets);
					}
					m_nOptNets = nCatnets;
					m_pOptNets = (CATNET<char, MAX_NODE_NAME, double>**) CATNET_MALLOC(
									nCatnets * sizeof(CATNET<char, MAX_NODE_NAME, double>*));

					if (m_pOptNets && pCatnets)
						memcpy(m_pOptNets, pCatnets, nCatnets * sizeof(CATNET<char,
							MAX_NODE_NAME, double>*));
					if (pCatnets)
						memset(pCatnets, 0, nCatnets * sizeof(CATNET<char,
							MAX_NODE_NAME, double>*));
					if (m_pOptOrder && m_pTestOrder[n])
						memcpy(m_pOptOrder, m_pTestOrder[n], m_numNodes * sizeof(int));
					stepaccept++;
					ordProb = newOrdProb[n];
				} else { 
					deltaLogLik = fLogLik - fOptLogLik;
					deltaLogLik += (newOrdProb[n] - ordProb);
					acceptprob = 0;
					if(tempCur > 0 && deltaLogLik <= 0)
						acceptprob = (double) exp((double)deltaLogLik/(double)tempCur);
					GetRNGstate();
					ftemp = (double) unif_rand();
					PutRNGstate();
					if (deltaLogLik > 0 || ftemp < acceptprob) {
						pOptNet = pCurNet;
						fOptLogLik = fLogLik;
						if (m_pOptNets) {
						    for (i = 0; i < m_nOptNets; i++) {
						        if (m_pOptNets[i]) {
						            delete m_pOptNets[i];
						        }
						    }
						    CATNET_FREE( m_pOptNets);
						}

						m_nOptNets = nCatnets;
						m_pOptNets = (CATNET<char, MAX_NODE_NAME, double>**) CATNET_MALLOC(
									nCatnets * sizeof(CATNET<char, MAX_NODE_NAME, double>*));
						if (m_pOptNets && pCatnets)
							memcpy(m_pOptNets, pCatnets, nCatnets * sizeof(CATNET<
								char, MAX_NODE_NAME, double>*));
						if (pCatnets)
							memset(pCatnets, 0, nCatnets * 
								sizeof(CATNET<char,MAX_NODE_NAME, double>*));
						if (m_pOptOrder && m_pTestOrder[n])
							memcpy(m_pOptOrder, m_pTestOrder[n], m_numNodes * sizeof(int));

						if (echo) {
							Rprintf("Accept order: ");
							if(m_numNodes < 32) {
								for (i = 0; i < m_numNodes; i++)
									Rprintf("%d ", m_pOptOrder[i]);
							}
							if(deltaLogLik > 0)
								Rprintf("\n");
							else
								Rprintf(" with prob %f < %f\n", ftemp, acceptprob);
						}
						ordProb = newOrdProb[n];

						stepaccept++;
						stepiters = n + 1;
						nstopCounter = 0;
					}
					else {
						deltaLogLik = 0;
					}

					if (deltaLogLik < 0)
						deltaLogLik = -deltaLogLik;
					if (deltaLogLik < stopDiff)
						nstopCounter++;
					if(nstopCounter >= tempCheckOrders)
						nstop = 1;
				}

			} // if(threadres == 0)

		} // for(n = 0; n < m_nDrives; n++)

		if(stepaccept > 0) {
			naccept++;
		}
		niter += stepiters;
		if (niter - nLastChangedTemp >= tempCheckOrders) {
			tempCur *= tempCoolFact;
			nLastChangedTemp = niter;
			if (echo) {
				Rprintf("Set temp to %f\n", tempCur);
			}
		}
		if (echo)
			Rprintf("Iteration %d\\%d\n", niter, maxIter);
					
		if (nstop)
			break;
	} //while(niter < maxIter)
	
	UNPROTECT(1); // rSamples

	if (!isNull(rParentSizes) && length(rParentSizes) == m_numNodes)
		UNPROTECT(1);
	if (!isNull(rPerturbations))
		UNPROTECT(1);
	if(!isNull(rNodeCats))
		UNPROTECT(1);
	if (!isNull(rParentsPool) && length(rParentsPool) == m_numNodes)
		UNPROTECT(1);
	if (!isNull(rFixedParentsPool) && length(rFixedParentsPool) == m_numNodes)
		UNPROTECT(1);
	if (!isNull(rMatEdgeLiks) && length(rMatEdgeLiks) == m_numNodes*m_numNodes)
		UNPROTECT(1);
	if (!isNull(rDirProbs) && length(rDirProbs) == m_numNodes*m_numNodes)
		UNPROTECT(1);

	if (echo && niter > 0) {
		Rprintf("Accepted\\Total Orders %d\\%d, (%f)\n", naccept, niter,
				(double) naccept / (double) niter);
	}

	if(newOrdProb)
		CATNET_FREE(newOrdProb);

	if (m_pTestOrder) {
		for (n = 0; n < m_nDrives; n++) {
			if (m_pTestOrder[n])
				CATNET_FREE( m_pTestOrder[n]);
		}
		CATNET_FREE( m_pTestOrder);
		m_pTestOrder = 0;
	}

	if (m_pTestOrderInverse) {
		for (n = 0; n < m_nDrives; n++) {
			if (m_pTestOrderInverse[n])
				CATNET_FREE( m_pTestOrderInverse[n]);
		}
		CATNET_FREE( m_pTestOrderInverse);
		m_pTestOrderInverse = 0;
	}

	for (n = 0; n < m_nDrives; n++) {
		if (m_pSearchParams && m_pSearchParams[n])
			delete m_pSearchParams[n];
		if (m_pDrives && m_pDrives[n])
			delete m_pDrives[n];
	}
	
	if(pNodeCondLiks) {
		CATNET_FREE(pNodeCondLiks);
		pNodeCondLiks = 0;
	}

	if (m_pSearchParams)
		CATNET_FREE(m_pSearchParams);
	m_pSearchParams = 0;
	if (m_pDrives)
		CATNET_FREE( m_pDrives);
	m_pDrives = 0;
	m_nDrives = 0;

	if (m_bUseCache) {
		MUTEX_DESTROY(m_cache_mutex);
	}
	
	if (!m_nOptNets || !m_pOptNets) {
		warning("No networks are found");
		return R_NilValue;
	}

	// create a R-list of catNetworks
	numnets = 0;
	for (i = 0; i < m_nOptNets; i++) {
		if (m_pOptNets[i]) {
			// set node names according to the right order
			//m_pOptNets[i]->setNodesOrder(m_pOptOrder);
			m_pOptNets[i]->setNodeNames(pNodeNames, m_pOptOrder);
			numnets++;
		}
	}

	if (m_pOptOrder)
		CATNET_FREE( m_pOptOrder);
	m_pOptOrder = 0;

	PROTECT(cnetlist = allocVector(VECSXP, numnets));

	inet = 0;
	for (i = 0; i < m_nOptNets; i++) {
		if (!m_pOptNets[i])
			continue;
		rcatnet = *m_pOptNets[i];
		PROTECT(cnetnode = rcatnet.genRcatnet("catNetwork"));
		SET_VECTOR_ELT(cnetlist, inet, cnetnode);
		UNPROTECT(1);
		inet++;
	}

	UNPROTECT(1); /* cnetlist */

	if (pNodeNames) {
		for (i = 0; i < m_numNodes; i++) {
			if (pNodeNames[i])
				CATNET_FREE(pNodeNames[i]);
			pNodeNames[i] = 0;
		}
		CATNET_FREE(pNodeNames);
		pNodeNames = 0;
	}

	return cnetlist;
}

