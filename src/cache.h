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
 * cache.h
 *
 *  Created on: Sep 21, 2009
 *      Author: nbalov
 */

/* 
 * version 1.15.1  12dec2016
 */

#ifndef CACHE_H
#define CACHE_H

#include "utils.h"
#include "problist.h"

template<class t_prob>
struct CATNET_CACHE_EL {
	int		nnode;
	PROB_LIST<t_prob> *pNodeProb;
	int		npars, *pPars;
	int		nPool, *pPool;
	t_prob		fLogLik;

	CATNET_CACHE_EL(int *ppool, int npool, int node, int *parset, int parsetsize, PROB_LIST<t_prob> *probNode, t_prob flik) {
		nPool = npool;
		pPool = (int*)CATNET_MALLOC(npool*sizeof(int));
		if (pPool)
			memcpy(pPool, ppool, npool*sizeof(int));
		nnode = node;
		npars = parsetsize;
		pPars = (int*)CATNET_MALLOC(npars*sizeof(int));
		if (pPars)
			memcpy(pPars, parset, npars*sizeof(int));
		pNodeProb = new PROB_LIST<t_prob>;
		if (pNodeProb && probNode)
			*pNodeProb = *probNode;
		fLogLik = flik;
	}

	~CATNET_CACHE_EL() {
		if(pPool)
			CATNET_FREE(pPool);
		if(pPars)
			CATNET_FREE(pPars);
		if(pNodeProb)
			delete pNodeProb;
	}
};

extern CATNET_CACHE_EL<double>		**g_pcache;
extern unsigned int			g_ncache;

void ReleaseCache();
void InitializeCache(int nslots, int ncachebits);

class c_cache {
protected:
	int m_numNodes, m_maxParentSet;
	int *m_pRorder, *m_pRorderInverse, *m_parBuff1, *m_parBuff2;
	int m_bUseCache, m_nCacheBits;

	void _release();
public:
	c_cache();
	~c_cache();

	void setCacheParams(int numNodes, int maxParentSet, int *pRorder, int *pRorderInverse);

	int getCachedProb(int *ppool, int poolsize, int node, int *parset, int parsize, 
						PROB_LIST<double> *probNode, double *pflik);
	int setCachedProb(int *ppool, int poolsize, int node, int *parset, int parsize, 
						PROB_LIST<double> *probNode, double flik);

};

#endif /* CACHE_H */
