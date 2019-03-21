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
#include "rcatnet_hist.h"
#include "rcatnet.h"

using namespace std;

THREAD_PROC_DEFINE(CatnetSearchHistThreadProc, pParam) {
	SEARCH_PARAMETERS *pSearchParams = (SEARCH_PARAMETERS*)pParam;
	if(pSearchParams && pSearchParams->m_pCaller) {
		CATNET_SEARCH2<char, MAX_NODE_NAME, double> *pCaller = 
			(CATNET_SEARCH2<char, MAX_NODE_NAME, double>*)pSearchParams->m_pCaller;
		pCaller->estimate(pSearchParams);
		pCaller->_exit_thread((void*)0);
	}
	return(0);
}

RCatnetSearchHist::RCatnetSearchHist() {
	m_nDrives = 0;
	m_pDrives = NULL;
	m_pTestOrder = NULL;
	m_pTestOrderInverse = NULL;
	m_bUseCache = 1;
	m_pSearchParams = NULL;
}

RCatnetSearchHist::~RCatnetSearchHist() {
	_release();
}

void RCatnetSearchHist::_release() {
	int i;

	if(m_pTestOrder) {
		for(i = 0; i < m_nDrives; i++) {
			if(m_pTestOrder[i])
				CATNET_FREE(m_pTestOrder[i]);
		}
		CATNET_FREE(m_pTestOrder);
		m_pTestOrder = NULL;
	}

	if(m_pTestOrderInverse) {
		for(i = 0; i < m_nDrives; i++) {
			if(m_pTestOrderInverse[i])
				CATNET_FREE(m_pTestOrderInverse[i]);
		}
		CATNET_FREE(m_pTestOrderInverse);
		m_pTestOrderInverse = NULL;
	}

	if(m_pSearchParams) {
		for(i = 0; i < m_nDrives; i++)
			if(m_pSearchParams[i])
				delete m_pSearchParams[i];
		CATNET_FREE(m_pSearchParams);
		m_pSearchParams = NULL;
	}

	if(m_pDrives) {
		for(i = 0; i < m_nDrives; i++) 
			if(m_pDrives[i])
				delete m_pDrives[i];
		CATNET_FREE(m_pDrives);
		m_pDrives = NULL;
		
	}
	m_nDrives = 0;
}

int *RCatnetSearchHist::_genOrder(const int *porder, int norder, int shuffles, int bjump) {
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

SEXP RCatnetSearchHist::search(SEXP rSamples, SEXP rPerturbations, 
			SEXP rMaxParents, SEXP rParentSizes, SEXP rMaxComplexity, SEXP rNodeCats, 
			SEXP rParentsPool, SEXP rFixedParentsPool, 
			SEXP rScore, SEXP rWeight, SEXP rMaxIter, 
			SEXP rThreads, SEXP rUseCache, SEXP rEcho) {

	int n, nn, i, j, k, len, maxComplexity, echo;
 	int *pRsamples, *pRperturbations, *pSamples, *pPerturbations, 
		**parentsPool, **fixedParentsPool, *pPool, *pParentSizes;
	double complx, fLogLik, ftemp, *pMatHisto;
	int niter, maxIter, nScoreSel, nWeight, nstop;

	MUTEX m_cache_mutex;

	int nCatnets, nCurNet;
	CATNET<char, MAX_NODE_NAME, double> **pCatnets, *pCurNet;
	const int **pParents, *pNumParents;

	SEXP dim, rnodecat, rparpool, rmat;

	_release();

	if(!isMatrix(rSamples))
		error("Data is not a matrix");

	PROTECT(rMaxParents = AS_INTEGER(rMaxParents));
	PROTECT(rMaxComplexity = AS_INTEGER(rMaxComplexity));
	PROTECT(rMaxIter = AS_INTEGER(rMaxIter));
	PROTECT(rScore = AS_CHARACTER(rScore));
	PROTECT(rWeight = AS_INTEGER(rWeight));
	PROTECT(rThreads = AS_INTEGER(rThreads));
	PROTECT(rUseCache = AS_LOGICAL(rUseCache));
	PROTECT(rEcho = AS_LOGICAL(rEcho));

	m_maxParentSet = INTEGER_POINTER(rMaxParents)[0];
	maxComplexity = INTEGER_POINTER(rMaxComplexity)[0];
	maxIter = INTEGER_POINTER(rMaxIter)[0];

	m_nDrives = INTEGER_POINTER(rThreads)[0];
	if(m_nDrives < 1)
		m_nDrives = 1;

	m_bUseCache = LOGICAL(rUseCache)[0];
	if(m_bUseCache && m_nDrives > 8) {
		// with many threads in parallel the cache is presumably inefficient - disable it!
		m_bUseCache = 0;
		warning("Cache is disabled, too many threads");
	}

	echo = LOGICAL(rEcho)[0];

	nScoreSel = 0;
	if(!strncmp(CHARACTER_VALUE(rScore), "AIC", 3))
		nScoreSel = 1;
	if(!strncmp(CHARACTER_VALUE(rScore), "BIC", 3))
		nScoreSel = 2;

	nWeight = INTEGER_POINTER(rWeight)[0];

	UNPROTECT(8);

	PROTECT(rSamples = AS_INTEGER(rSamples));
	pRsamples = INTEGER(rSamples);

	dim = GET_DIM(rSamples);
	m_numNodes = INTEGER(dim)[0];
	m_numSamples = INTEGER(dim)[1];

	// create a R histogram matrix
	PROTECT(rmat = NEW_NUMERIC(m_numNodes*m_numNodes));
	pMatHisto = NUMERIC_POINTER(rmat);
	if(!pMatHisto)
		return R_NilValue;
	memset(pMatHisto, 0, m_numNodes*m_numNodes*sizeof(double));

	m_pTestOrder = (int**)CATNET_MALLOC(m_nDrives*sizeof(int*));
	if (m_pTestOrder)
		memset(m_pTestOrder, 0, m_nDrives*sizeof(int*));
	m_pTestOrderInverse = (int**)CATNET_MALLOC(m_nDrives*sizeof(int*));
	for(n = 0; n < m_nDrives; n++) {
		m_pTestOrderInverse[n] = (int*)CATNET_MALLOC(m_numNodes*sizeof(int));
	}

	m_pDrives = (CATNET_SEARCH2<char, MAX_NODE_NAME, double>**)CATNET_MALLOC(
			m_nDrives*sizeof(CATNET_SEARCH2<char, MAX_NODE_NAME, double>*));
	for(n = 0; n < m_nDrives; n++) {
		m_pDrives[n] = new CATNET_SEARCH2<char, MAX_NODE_NAME, double>();
	}

	if(m_bUseCache) {
		MUTEX_INIT(m_cache_mutex);
	}

	if(!isNull(rParentSizes) && length(rParentSizes) == m_numNodes)
		PROTECT(rParentSizes = AS_INTEGER(rParentSizes));
	if(!isNull(rPerturbations))
		PROTECT(rPerturbations = AS_INTEGER(rPerturbations));
	if(!isNull(rNodeCats))
		PROTECT(rNodeCats = AS_LIST(rNodeCats));
	if(!isNull(rParentsPool) && length(rParentsPool) == m_numNodes)
		PROTECT(rParentsPool = AS_LIST(rParentsPool));
	if(!isNull(rFixedParentsPool) && length(rFixedParentsPool) == m_numNodes)
		PROTECT(rFixedParentsPool = AS_LIST(rFixedParentsPool));

	m_pSearchParams = (SEARCH_PARAMETERS**)CATNET_MALLOC(m_nDrives*sizeof(SEARCH_PARAMETERS*));
	for(n = 0; n < m_nDrives; n++) {
		m_pSearchParams[n] = new SEARCH_PARAMETERS (
			m_numNodes, m_numSamples, m_maxParentSet, maxComplexity, 
			(m_nDrives < 2)?echo:0, 
			!isNull(rNodeCats), 
			!isNull(rParentSizes), !isNull(rPerturbations), 
			!isNull(rParentsPool), !isNull(rFixedParentsPool), FALSE, 0, 
			m_bUseCache?&m_cache_mutex:0, m_pDrives[n]);
	}

	nstop = 0;

	niter = 0;
	while(niter < maxIter) {

		for(n = 0; n < m_nDrives; n++)  {

		if(m_pTestOrder[n])
			CATNET_FREE(m_pTestOrder[n]);
		m_pTestOrder[n] = _genOrder(0, m_numNodes, 0, 0);
		if(!m_pTestOrder[n])
			error("_genOrder returns an error");

		nstop = 0;
		if(n > 0) {
			for(nn = 0; nn < n; nn++) {
				nstop = 1;
				for(i = 0; i < m_numNodes; i++) 
					if(m_pTestOrder[nn][i] != m_pTestOrder[n][i]) {
						nstop = 0;
						break;
					}
				if(nstop)
					break;
			}
		}
		if(nstop) {
			Rprintf("repeated order in %d\n", n);
			continue;
		}

		for(i = 0; i < m_numNodes; i++)
			m_pTestOrderInverse[n][m_pTestOrder[n][i]-1] = i + 1;

		if(m_bUseCache)
			m_pDrives[n]->setCacheParams(m_numNodes, m_maxParentSet, 
							m_pTestOrder[n], m_pTestOrderInverse[n]);

		if(!isNull(rParentSizes) && length(rParentSizes) == m_numNodes) {
			pParentSizes = m_pSearchParams[n]->m_pParentSizes;
			if (pParentSizes && INTEGER(rParentSizes))
				memcpy(pParentSizes, INTEGER(rParentSizes), m_numNodes*sizeof(int));
			for(i = 0; i < m_numNodes; i++)
				pParentSizes[i] = INTEGER(rParentSizes)[m_pTestOrder[n][i] - 1];
		}

		pSamples = m_pSearchParams[n]->m_pSamples;
		for(j = 0; j < m_numSamples; j++) {
			for(i = 0; i < m_numNodes; i++) {
				pSamples[j*m_numNodes + i] = pRsamples[j*m_numNodes + m_pTestOrder[n][i] - 1];
				if(R_IsNA(pSamples[j*m_numNodes + i]) || pSamples[j*m_numNodes + i] < 1)
					pSamples[j*m_numNodes + i] = CATNET_NAN;
			}
		}

		if(!isNull(rPerturbations)) {
			pPerturbations = m_pSearchParams[n]->m_pPerturbations;
			pRperturbations = INTEGER(rPerturbations);
			for(j = 0; j < m_numSamples; j++) {
				for(i = 0; i < m_numNodes; i++)
					pPerturbations[j*m_numNodes + i] = pRperturbations[j*m_numNodes + m_pTestOrder[n][i] - 1];
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

		if(!isNull(rParentsPool) && length(rParentsPool) == m_numNodes) {
			parentsPool = m_pSearchParams[n]->m_parentsPool;
			for(i = 0; i < m_numNodes; i++) {
				rparpool = AS_INTEGER(VECTOR_ELT(rParentsPool, (int)(m_pTestOrder[n][i] - 1)));
				len = length(rparpool);
				if(isVector(rparpool) && len > 0 && len <= m_numNodes) {
					pPool = INTEGER(rparpool);
					for(j = 0; j < len; j++) {
						if(pPool[j] > 0 && pPool[j] <= m_numNodes) {
							for(k = 0; k < m_numNodes; k++)
								if(pPool[j] == m_pTestOrder[n][k])
									break;
							if(k < m_numNodes)
								parentsPool[i][j] = k;
							else
								parentsPool[i][j] = -1;
						}
						else
							parentsPool[i][j] = -1;
					}
					for(; j < m_numNodes; j++)
						parentsPool[i][j] = -1;
				}
				else {
					for(j = 0; j < m_numNodes; j++)
						parentsPool[i][j] = -1;
				}
			}
		}

		if(!isNull(rFixedParentsPool) && length(rFixedParentsPool) == m_numNodes) {
			fixedParentsPool = m_pSearchParams[n]->m_fixedParentsPool;
			for(i = 0; i < m_numNodes; i++) {
				rparpool = AS_INTEGER(VECTOR_ELT(rFixedParentsPool, (int)(m_pTestOrder[n][i] - 1)));
				len = length(rparpool);
				if(isVector(rparpool) && len > 0 && len <= m_numNodes) {
				 	if(m_maxParentSet < len)
				    		m_maxParentSet = len;
					pPool = INTEGER(rparpool);
					for(j = 0; j < len; j++) {
						if(pPool[j] > 0 && pPool[j] <= m_numNodes) {
							for(k = 0; k < m_numNodes; k++)
								if(pPool[j] == m_pTestOrder[n][k])
									break;
							if(k < m_numNodes)
								fixedParentsPool[i][j] = k;
							else
								fixedParentsPool[i][j] = -1;
						}
						else
							fixedParentsPool[i][j] = -1;
					}
					for(; j < m_numNodes; j++)
						fixedParentsPool[i][j] = -1;
				}
				else {
					for(j = 0; j < m_numNodes; j++)
						fixedParentsPool[i][j] = -1;
				}
			}
		}

		m_pDrives[n] -> _start_thread(CatnetSearchHistThreadProc, (void*)m_pSearchParams[n]);

		} // for(n = 0; n < m_nDrives; n++)

		nstop = 1;

		for(n = 0; n < m_nDrives; n++)  {

		if(!m_pDrives[n] ->_is_running())
			continue;

		nstop = 0;

		int threadres = m_pDrives[n] -> _join_thread();
		threadres = m_pDrives[n] -> _stop_thread();

		nCatnets = 0;
		if(threadres == 0) {
			nCatnets = m_pDrives[n]->numCatnets();
			pCatnets = m_pDrives[n]->catnets();
		}
		if(nCatnets > 0 && pCatnets) {

			fLogLik = 0;
			pCurNet = NULL;
			switch(nScoreSel) {
			case 1: // AIC
				fLogLik = -FLT_MAX;
				for(nCurNet = 0; nCurNet < nCatnets; nCurNet++) {
					if(!pCatnets[nCurNet])
						continue;
					// net complexity = nCurNet + m_numNodes*(1 + maxCategories)
					complx = m_numSamples*pCatnets[nCurNet]->loglik() - pCatnets[nCurNet]->complexity();
					if(complx > fLogLik) {
						fLogLik = complx;
						pCurNet = pCatnets[nCurNet];
					}
				}
			break;
			case 2:  // BIC
				fLogLik = -FLT_MAX;
				ftemp = 0.5*log((double)m_numSamples);
				for(nCurNet = 0; nCurNet < nCatnets; nCurNet++) {
					if(!pCatnets[nCurNet])
						continue;
					complx = m_numSamples*pCatnets[nCurNet]->loglik() - ftemp * pCatnets[nCurNet]->complexity();
					if(complx > fLogLik) {
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
			if(!pCurNet) {
				continue;
			}

			if(nWeight <= 0)
				ftemp = 1;
			else if(nWeight == 1)
				ftemp =(double)exp( pCurNet->loglik() / (double)(m_numSamples*m_numNodes));
			else
				ftemp = (double)exp((double)fLogLik / (double)(m_numSamples*m_numNodes));

			pNumParents = pCurNet->numParents();
			pParents = pCurNet->parents();
			for(j = 0; j < m_numNodes; j++) {
				if(pNumParents[j] <= 0)
					continue;
				for(i = 0; i < pNumParents[j]; i++) {
				// the matrix rmat is in [parent, child] format, 
				// which in R is visible in [child, parent] matrix format
					pMatHisto[
						(m_pTestOrder[n][pParents[j][i]]-1)*m_numNodes + 
						 m_pTestOrder[n][j] - 1] += ftemp;
				}
			}

		} // if(threadres == 0)

		niter++;
		if(echo) {
			Rprintf("Iteration %d\\%d\n", niter, maxIter);
		}
		
		} // for(n = 0; n < m_nDrives; n++)

		if(nstop)
			break;
	} //while(niter < maxIter)

	UNPROTECT(1); // rSamples
	if(!isNull(rParentSizes) && length(rParentSizes) == m_numNodes) 
		UNPROTECT(1);
	if(!isNull(rPerturbations))
		UNPROTECT(1);
	if(!isNull(rNodeCats))
		UNPROTECT(1);
	if(!isNull(rParentsPool) && length(rParentsPool) == m_numNodes)
		UNPROTECT(1);
	if(!isNull(rFixedParentsPool) && length(rFixedParentsPool) == m_numNodes)
		UNPROTECT(1);

	if(m_pTestOrder) {
		for(n = 0; n < m_nDrives; n++) {
			if(m_pTestOrder[n])
				CATNET_FREE(m_pTestOrder[n]);
		}
		CATNET_FREE(m_pTestOrder);
		m_pTestOrder = 0;
	}

	if(m_pTestOrderInverse) {
		for(n = 0; n < m_nDrives; n++) {
			if(m_pTestOrderInverse[n])
				CATNET_FREE(m_pTestOrderInverse[n]);
		}
		CATNET_FREE(m_pTestOrderInverse);
		m_pTestOrderInverse = 0;
	}

	for(n = 0; n < m_nDrives; n++) {
		if(m_pSearchParams && m_pSearchParams[n])
			delete m_pSearchParams[n];
		if(m_pDrives && m_pDrives[n])
			delete m_pDrives[n];
	}

	if(m_pSearchParams)
		CATNET_FREE(m_pSearchParams);
	m_pSearchParams = 0;
	
	if(m_pDrives)
		CATNET_FREE(m_pDrives);
	m_pDrives = 0;
	m_nDrives = 0;

	if(m_bUseCache) {
		MUTEX_DESTROY(m_cache_mutex);
	}

	UNPROTECT(1);//rmat

	return rmat;
}

