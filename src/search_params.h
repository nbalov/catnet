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
 * search_params.h
 *
 *  Created on: Sep 25, 2009
 *      Author: nbalov
 */

/* 
 * version 1.15.1  12dec2016
 */

#include "utils.h"
#include "catnet_class.h"
#include "thread.h"
#include "cache.h"

#ifndef SEARCH_PARAMS_H
#define SEARCH_PARAMS_H

struct SEARCH_PARAMETERS {
	int m_numNodes;
	int m_numSamples;
	int *m_pSamples;
	int *m_pNodeNumCats;
	int **m_pNodeCats;
	int *m_pPerturbations;
	int m_maxParentSet;
	int *m_pParentSizes;
	int m_maxComplexity;
	int m_maxParentsPool;
	int **m_parentsPool;
	int **m_fixedParentsPool;
	double *m_matEdgeLiks;
	double *m_matNodeCondLiks;
	int m_echo;
	MUTEX *m_pCacheMutex;
	int m_seed;

	void *m_pCaller;

	SEARCH_PARAMETERS(int numNodes, int numSamples, int maxParentSet, int maxComplexity, int echo, 
			int hasCats = 0, int hasParentSizes = 0, int hasPerturbations = 0, 
			int hasParentsPool = 0, int hasFixedParentsPool = 0, 
			int hasEdgeLiks = 0, int maxParentsPool=0, 
			MUTEX *pCacheMutex = 0, void *pCaller = 0, int nSeed = 0) {
		int i;
		m_numNodes = numNodes;
		m_numSamples = numSamples;
		
		m_maxParentSet = maxParentSet;
		m_maxComplexity = maxComplexity;
		
		m_echo = echo;
		m_pCacheMutex = pCacheMutex;
		m_pCaller = pCaller;

		m_pParentSizes = 0;
		if(hasParentSizes) { 
			m_pParentSizes = (int*)CATNET_MALLOC(m_numNodes*sizeof(int));
			if (m_pParentSizes) {
				for(i = 0; i < m_numNodes; i++) 
					m_pParentSizes[i] = m_maxParentSet;
			}
		}

		m_pSamples = (int*)CATNET_MALLOC(m_numNodes*m_numSamples*sizeof(int));

		m_pNodeNumCats = 0;
		m_pNodeCats = 0;
		if(hasCats) { 
			m_pNodeNumCats = (int*)CATNET_MALLOC(m_numNodes*sizeof(int));
			if (m_pNodeNumCats)
				memset(m_pNodeNumCats, 0, m_numNodes*sizeof(int));
			m_pNodeCats = (int**)CATNET_MALLOC(m_numNodes*sizeof(int*));
			if (m_pNodeCats)
				memset(m_pNodeCats, 0, m_numNodes*sizeof(int*));
		}

		m_pPerturbations = 0;
		if(hasPerturbations)
			m_pPerturbations = (int*)CATNET_MALLOC(m_numNodes*m_numSamples*sizeof(int));

		m_maxParentsPool = maxParentsPool;
		m_matNodeCondLiks = 0;
		/* apply the maxParentsPool restriction only if there are no user specified parent pools */
		if(!hasParentsPool && m_maxParentsPool >= 1) {
			m_matNodeCondLiks = (double*)CATNET_MALLOC(m_numNodes*m_numNodes*sizeof(double));
			hasParentsPool = 1;
		}

		m_parentsPool = 0;
		if(hasParentsPool) {
			m_parentsPool = (int**)CATNET_MALLOC(m_numNodes*sizeof(int*));
			if (m_parentsPool) {
				for(i = 0; i < m_numNodes; i++) {
					m_parentsPool[i] = (int*)CATNET_MALLOC(m_numNodes*sizeof(int));
				}
			}
		}

		m_fixedParentsPool = 0;
		if(hasFixedParentsPool) {
			m_fixedParentsPool = (int**)CATNET_MALLOC(m_numNodes*sizeof(int*));
			if (m_fixedParentsPool) {
				for(i = 0; i < m_numNodes; i++) {
					m_fixedParentsPool[i] = (int*)CATNET_MALLOC(m_numNodes*sizeof(int));
				}
			}
		}
		
		m_matEdgeLiks = 0;
		if(hasEdgeLiks) {
			m_matEdgeLiks = (double*)CATNET_MALLOC(m_numNodes*m_numNodes*sizeof(double));
		}

		m_seed = nSeed;
	}

	~SEARCH_PARAMETERS() {
		int i;
		if(m_pParentSizes)
			CATNET_FREE(m_pParentSizes);
		if(m_pSamples)
			CATNET_FREE(m_pSamples);
		if(m_pPerturbations)
			CATNET_FREE(m_pPerturbations);
		if(m_pNodeCats) {
			for(i = 0; i < m_numNodes; i++)
				if(m_pNodeCats[i])
					CATNET_FREE(m_pNodeCats[i]);
			CATNET_FREE(m_pNodeCats);
		}
		if(m_pNodeNumCats)
			CATNET_FREE(m_pNodeNumCats);
		if(m_parentsPool) {
			for(i = 0; i < m_numNodes; i++)
				if(m_parentsPool[i])
					CATNET_FREE(m_parentsPool[i]);
			CATNET_FREE(m_parentsPool);
		}
		if(m_fixedParentsPool) {
			for(i = 0; i < m_numNodes; i++)
				if(m_fixedParentsPool[i])
					CATNET_FREE(m_fixedParentsPool[i]);
			CATNET_FREE(m_fixedParentsPool);
		}
		if(m_matEdgeLiks) {
			CATNET_FREE(m_matEdgeLiks);
		}
		if(m_matNodeCondLiks) {
			CATNET_FREE(m_matNodeCondLiks);
		}
	}
};

#endif /* SEARCH_PARAMS_H */

