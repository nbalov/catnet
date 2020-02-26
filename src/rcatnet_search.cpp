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
 * rcatnet_search.cpp
 *
 *  Created on: Sep 21, 2009
 *      Author: nbalov
 */

/* 
 * version 1.15.1  12dec2016
 */

#include "utils.h"
#include "catnet_class.h"
#include "rcatnet_search.h"
#include "rcatnet.h"

using namespace std;

RCatnetSearch::RCatnetSearch() {
	m_pRorder = 0;
	m_pRorderInverse = 0;
	m_bUseCache = 1;
	m_pSearchParams = 0;
}

RCatnetSearch::~RCatnetSearch() {
	if(m_pRorder)
		CATNET_FREE(m_pRorder);
	m_pRorder = 0;
	if(m_pRorderInverse)
		CATNET_FREE(m_pRorderInverse);
	m_pRorderInverse = 0;
	if(m_pSearchParams)
		delete m_pSearchParams;
	m_pSearchParams = 0;
}

SEXP RCatnetSearch::estimateCatnets(SEXP rSamples, SEXP rPerturbations, 
			SEXP rMaxParents, SEXP rParentSizes, SEXP rMaxComplexity, SEXP rOrder, SEXP rNodeCats, 
			SEXP rParentsPool, SEXP rFixedParentsPool, SEXP rMatEdgeLiks, SEXP rUseCache, SEXP rEcho) {

	int i, j, k, len, maxComplexity, numnets, inet, echo;
 	int *pRsamples, *pRperturbations, *pSamples, *pPerturbations, 
		**parentsPool, **fixedParentsPool, *pPool, *pParentSizes;
	double *matEdgeLiks, *pMatEdgeLiks;

	RCatnet rcatnet;
	SEXP dim, rnodecat, rparpool, cnetlist, cnetnode;

	if(!isMatrix(rSamples))
		error("Data is not a matrix");

	PROTECT(rMaxParents = AS_INTEGER(rMaxParents));
	m_maxParentSet = INTEGER_POINTER(rMaxParents)[0];
	UNPROTECT(1);

	PROTECT(rMaxComplexity = AS_INTEGER(rMaxComplexity));
	maxComplexity = INTEGER_POINTER(rMaxComplexity)[0];
	UNPROTECT(1);

	PROTECT(rUseCache = AS_LOGICAL(rUseCache));
	m_bUseCache = LOGICAL(rUseCache)[0];
	UNPROTECT(1);

	PROTECT(rEcho = AS_LOGICAL(rEcho));
	echo = LOGICAL(rEcho)[0];
	UNPROTECT(1);

	PROTECT(rSamples = AS_INTEGER(rSamples));
	dim = GET_DIM(rSamples);
	m_numNodes = INTEGER(dim)[0];
	m_numSamples = INTEGER(dim)[1]; 

	if(m_pRorder)
		CATNET_FREE(m_pRorder);
	m_pRorder = (int*)CATNET_MALLOC(m_numNodes*sizeof(int));
	if(m_pRorderInverse)
		CATNET_FREE(m_pRorderInverse);
	m_pRorderInverse = (int*)CATNET_MALLOC(m_numNodes*sizeof(int));

	if(!m_pRorder || !m_pRorderInverse) {
		if(m_pRorder)
			CATNET_FREE(m_pRorder);
		m_pRorder = 0;
		if(m_pRorderInverse)
			CATNET_FREE(m_pRorderInverse);
		m_pRorderInverse = 0;
		UNPROTECT(1);
		return R_NilValue;
	}

	PROTECT(rOrder = AS_INTEGER(rOrder));
	if(length(rOrder) < m_numNodes) {
		warning("Invalid nodeOrder parameter - reset to default node order.");
		for(i = 0; i < m_numNodes; i++)
			m_pRorder[i] = i + 1;
	}
	else {
		if (INTEGER(rOrder))
			memcpy(m_pRorder, INTEGER(rOrder), m_numNodes*sizeof(int));
	}
	for(i = 0; i < m_numNodes; i++) {
		if(m_pRorder[i] <= 0 || m_pRorder[i] > m_numNodes) {
			error("Invalid nodeOrder parameter");		
		}	
		else
			m_pRorderInverse[m_pRorder[i]-1] = i + 1;
	}
	UNPROTECT(1);

	if(m_pSearchParams)
		delete m_pSearchParams;
	m_pSearchParams = new SEARCH_PARAMETERS(
		m_numNodes, m_numSamples, 
		m_maxParentSet, maxComplexity, echo, 
		!isNull(rNodeCats), 
		!isNull(rParentSizes), !isNull(rPerturbations), 
		!isNull(rParentsPool), !isNull(rFixedParentsPool), 
		!isNull(rMatEdgeLiks), 0, 
		NULL, this);

	if(!isNull(rParentSizes)) {
		pParentSizes = m_pSearchParams->m_pParentSizes;
		PROTECT(rParentSizes = AS_INTEGER(rParentSizes));
		if(length(rParentSizes) == m_numNodes) { 
			for(i = 0; i < m_numNodes; i++)
				pParentSizes[i] = INTEGER(rParentSizes)[m_pRorder[i] - 1];
		}
		UNPROTECT(1);
	}

	pSamples = m_pSearchParams->m_pSamples;
	pRsamples = INTEGER(rSamples);
	for(j = 0; j < m_numSamples; j++) {
		for(i = 0; i < m_numNodes; i++) {
			pSamples[j*m_numNodes + i] = pRsamples[j*m_numNodes + m_pRorder[i] - 1];
			if(R_IsNA(pSamples[j*m_numNodes + i]) || pSamples[j*m_numNodes + i] < 1) {
				pSamples[j*m_numNodes + i] = CATNET_NAN;
			}
		}
	}
	UNPROTECT(1); // rSamples

	pPerturbations = 0;
	if(!isNull(rPerturbations)) {
		PROTECT(rPerturbations = AS_INTEGER(rPerturbations));
		pPerturbations = m_pSearchParams->m_pPerturbations;
		pRperturbations = INTEGER_POINTER(rPerturbations);
		for(j = 0; j < m_numSamples; j++) {
			for(i = 0; i < m_numNodes; i++) {
				pPerturbations[j*m_numNodes + i] = pRperturbations[j*m_numNodes + m_pRorder[i] - 1];
			}
		}
		UNPROTECT(1);
	}

	if(!isNull(rNodeCats)) {
		PROTECT(rNodeCats = AS_LIST(rNodeCats));
		for(i = 0; i < m_numNodes; i++) {
			rnodecat = AS_INTEGER(VECTOR_ELT(rNodeCats, (int)(m_pRorder[i] - 1)));
			len = length(rnodecat);
			if(isVector(rnodecat) && len > 0) {
				m_pSearchParams->m_pNodeNumCats[i] = len;
				m_pSearchParams->m_pNodeCats[i] = (int*)CATNET_MALLOC(len*sizeof(int));
				for(j = 0; j < len; j++)
					m_pSearchParams->m_pNodeCats[i][j] = INTEGER(rnodecat)[j];
			}
		}
		UNPROTECT(1);
	}

	parentsPool = 0;
	if(!isNull(rParentsPool)) {
		PROTECT(rParentsPool = AS_LIST(rParentsPool));
		parentsPool = m_pSearchParams->m_parentsPool;
		for(i = 0; i < m_numNodes; i++) {
			rparpool = AS_INTEGER(VECTOR_ELT(rParentsPool, (int)(m_pRorder[i] - 1)));
			len = length(rparpool);
			if(isVector(rparpool) && len > 0 && len <= m_numNodes) {
				pPool = INTEGER(rparpool);
				for(j = 0; j < len; j++) {
					if(pPool[j] > 0 && pPool[j] <= m_numNodes) {
						for(k = 0; k < m_numNodes; k++)
							if(pPool[j] == m_pRorder[k])
								break;
						if(k < m_numNodes)
							parentsPool[i][j] = k;
						else
							parentsPool[i][j] = -1;
					}
				}
				for(; j < m_numNodes; j++)
					parentsPool[i][j] = -1;
			}
			else {
				//error("Invalid parentsPool parameter");
				for(j = 0; j < m_numNodes; j++)
					parentsPool[i][j] = -1;
				// no parents allow - achieved by a reference to itself
				parentsPool[i][0] = i;
			}
		}
		UNPROTECT(1);
	}

	fixedParentsPool = 0;
	if(!isNull(rFixedParentsPool)) {
		PROTECT(rFixedParentsPool = AS_LIST(rFixedParentsPool));
		fixedParentsPool = m_pSearchParams->m_fixedParentsPool;
		for(i = 0; i < m_numNodes; i++) {
			rparpool = AS_INTEGER(VECTOR_ELT(rFixedParentsPool, (int)(m_pRorder[i] - 1)));
			len = length(rparpool);
			if(isVector(rparpool) && len > 0 && len <= m_numNodes) {
			 	if(m_maxParentSet < len)
			    		m_maxParentSet = len;
				pPool = INTEGER(rparpool);
				for(j = 0; j < len; j++) {
					if(pPool[j] > 0 && pPool[j] <= m_numNodes) {
						for(k = 0; k < m_numNodes; k++)
							if(pPool[j] == m_pRorder[k])
								break;
						if(k < m_numNodes)
							fixedParentsPool[i][j] = k;
						else
							fixedParentsPool[i][j] = -1;
					}
				}
				for(; j < m_numNodes; j++)
					fixedParentsPool[i][j] = -1;
			}
			else {
				for(j = 0; j < m_numNodes; j++)
					fixedParentsPool[i][j] = -1;
			}
		}
		UNPROTECT(1);
	}

	if(!isNull(rMatEdgeLiks) && m_pSearchParams->m_matEdgeLiks) {
		PROTECT(rMatEdgeLiks = AS_NUMERIC(rMatEdgeLiks));
		matEdgeLiks = m_pSearchParams->m_matEdgeLiks;
		pMatEdgeLiks = REAL(rMatEdgeLiks);
		for(j = 0; j < m_numNodes; j++) {
			for(i = 0; i < m_numNodes; i++) {
				matEdgeLiks[j*m_numNodes + i] = pMatEdgeLiks[(m_pRorder[j] - 1)*m_numNodes + m_pRorder[i] - 1];
			}
		}
		UNPROTECT(1);
	}

	if(m_bUseCache)
		setCacheParams(m_numNodes, m_maxParentSet, m_pRorder, m_pRorderInverse);

	estimate(m_pSearchParams);

	if(m_pSearchParams)
		delete m_pSearchParams;
	m_pSearchParams = 0;

	if(!m_nCatnets || !m_pCatnets) {
		warning("No networks are found");
		return R_NilValue;
	}

	// create a R-list of catNetworks
	numnets = 0;
	for(i = 0; i < m_nCatnets; i++) {
		if(m_pCatnets[i]) {
			m_pCatnets[i]->setNodesOrder(m_pRorder);
			numnets++;
		}
	}

	PROTECT(cnetlist = allocVector(VECSXP, numnets));
	inet = 0;
	for(i = 0; i < m_nCatnets; i++) {
		if(!m_pCatnets[i])
			continue;

		rcatnet = *m_pCatnets[i];

		//rcatnet.setCategoryIndices(m_pNodeNumCats, m_pNodeCats);
		PROTECT(cnetnode = rcatnet.genRcatnet("catNetwork"));

		SET_VECTOR_ELT(cnetlist, inet, cnetnode);
		UNPROTECT(1);
		inet++;
	}

	UNPROTECT(1);//cnetlist

	return cnetlist;
}
