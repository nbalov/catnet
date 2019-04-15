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
 * cache.cpp
 *
 *  Created on: Sep 21, 2009
 *      Author: nbalov
 */

/* 
 * version 1.15.1  12dec2016
 */

#include "cache.h"

#define CACHE_MAX_LEN_POOL 50
int PRIMES_NUM = 180;
unsigned PRIMES_1000[] = {
      2,      3,      5,      7,     11,     13,     17,     19,     23,     29, 
     31,     37,     41,     43,     47,     53,     59,     61,     67,     71, 
     73,     79,     83,     89,     97,    101,    103,    107,    109,    113, 
    127,    131,    137,    139,    149,    151,    157,    163,    167,    173, 
    179,    181,    191,    193,    197,    199,    211,    223,    227,    229, 
    233,    239,    241,    251,    257,    263,    269,    271,    277,    281, 
    283,    293,    307,    311,    313,    317,    331,    337,    347,    349, 
    353,    359,    367,    373,    379,    383,    389,    397,    401,    409, 
    419,    421,    431,    433,    439,    443,    449,    457,    461,    463, 
    467,    479,    487,    491,    499,    503,    509,    521,    523,    541, 
    547,    557,    563,    569,    571,    577,    587,    593,    599,    601, 
    607,    613,    617,    619,    631,    641,    643,    647,    653,    659, 
    661,    673,    677,    683,    691,    701,    709,    719,    727,    733, 
    739,    743,    751,    757,    761,    769,    773,    787,    797,    809, 
    811,    821,    823,    827,    829,    839,    853,    857,    859,    863, 
    877,    881,    883,    887,    907,    911,    919,    929,    937,    941, 
    947,    953,    967,    971,    977,    983,    991,    997,   1009,   1013, 
   1019,   1021,   1031,   1033,   1039,   1049,   1051,   1061,   1063,   1069 
};

//#define DEBUG_INFO

CATNET_CACHE_EL<double>		**g_pcache = 0;
unsigned int			g_ncache = 0;
unsigned int			g_nCacheBits = 0;

void ReleaseCache() {
	unsigned i;
	if(g_pcache && g_ncache > 0) {
		for(i = 0; i < g_ncache; i++) {
			if(g_pcache[i])
				delete g_pcache[i];
			g_pcache[i] = 0;
		}
	}
	if(g_pcache)
		CATNET_FREE(g_pcache);
	g_pcache = 0;
	g_ncache = 0;
	g_nCacheBits = 0;
}
	
void InitializeCache(int nslots, int ncachebits) {
	ReleaseCache();
	if(nslots < 1)
		nslots = 1;
	g_ncache = nslots;
	g_pcache = (CATNET_CACHE_EL<double>**)CATNET_MALLOC(g_ncache*sizeof(CATNET_CACHE_EL<double>*));
	if (g_pcache)
		memset(g_pcache, 0, g_ncache*sizeof(CATNET_CACHE_EL<double>*));
	g_nCacheBits = ncachebits;
}


c_cache::c_cache() {
	m_numNodes = 0;
	m_pRorder = 0;
	m_pRorderInverse = 0;
	m_parBuff1 = 0;
	m_parBuff2 = 0;
	m_bUseCache = 0;
}

void c_cache::_release() {
	if(m_pRorder)
		CATNET_FREE(m_pRorder);
	m_pRorder = 0;
	if(m_pRorderInverse)
		CATNET_FREE(m_pRorderInverse);
	m_pRorderInverse = 0;
	if(m_parBuff1)
		CATNET_FREE(m_parBuff1);
	m_parBuff1 = 0;
	if(m_parBuff2)
		CATNET_FREE(m_parBuff2);
	m_parBuff2 = 0;
}

c_cache::~c_cache() {
	_release();
}


void c_cache::setCacheParams(int numNodes, int maxParentSet, int *pRorder, int *pRorderInverse) {
	if(numNodes < 1 || maxParentSet < 0 || !pRorder || !pRorderInverse)
		return;
	if(m_numNodes != numNodes) {
		_release();
	}
	m_numNodes = numNodes;
	m_maxParentSet = maxParentSet;
	if(!m_pRorder)
		m_pRorder = (int*)CATNET_MALLOC(m_numNodes*sizeof(int));
	if(m_pRorder && pRorder)
		memcpy(m_pRorder, pRorder, m_numNodes*sizeof(int));
	if(!m_pRorderInverse)
		m_pRorderInverse = (int*)CATNET_MALLOC(m_numNodes*sizeof(int));
	if (m_pRorderInverse && pRorderInverse)
		memcpy(m_pRorderInverse, pRorderInverse, m_numNodes*sizeof(int));
	if(!m_parBuff1)
		m_parBuff1 = (int*)CATNET_MALLOC(m_numNodes*sizeof(int));
	if(!m_parBuff2)
		m_parBuff2 = (int*)CATNET_MALLOC(m_numNodes*sizeof(int));
	m_bUseCache = 1;
}

int c_cache::getCachedProb(int *ppool, int poolsize, int node, int *parset, int parsize, 
				PROB_LIST<double> *probNode, double *pflik) {

		int i, j;
		unsigned nlookup;
		CATNET_CACHE_EL<double> *pCacheEl = 0;

		if(!m_bUseCache)
			return 0;
		if(m_numNodes < 1 || !m_pRorder)
			return 0;
		if(!g_pcache)
			return 0;
#ifdef DEBUG_INFO
char str[256];
sprintf(str,"getCachedProb node=%d, pool = ", node);Rprintf(str);
for(i = 0; i < poolsize; i++) {
	sprintf(str,"%d  ", ppool[i]);Rprintf(str);
}Rprintf("\n");
#endif
		node = m_pRorder[node];
		for(i = 0; i < poolsize; i++)
			m_parBuff1[i] = m_pRorder[ppool[i]];
		_quick_sort<int>(m_parBuff1, poolsize);

#ifdef DEBUG_INFO
sprintf(str,"    reordered node=%d, pool = ", node);Rprintf(str);
for(i = 0; i < poolsize; i++) {
	sprintf(str,"%d  ", m_parBuff1[i]); Rprintf(str);
}Rprintf("\n");
#endif
		nlookup = 1;
		for(i = 0; i < poolsize; i++) {
			j = m_parBuff1[i] - 1;
			while(j >= PRIMES_NUM)
				j -= PRIMES_NUM;
			nlookup *= PRIMES_1000[PRIMES_NUM - j - 1];
			while(nlookup >= g_ncache)
				nlookup -= g_ncache;		
		}
		nlookup = (nlookup << g_nCacheBits) + node + m_numNodes*parsize;
		while(nlookup >= g_ncache)
			nlookup -= g_ncache;

#ifdef DEBUG_INFO
sprintf(str,"nlookup=%d\n", nlookup);Rprintf(str);
#endif

		pCacheEl = g_pcache[nlookup];
		if(!pCacheEl) 
			return 0;
		if(pCacheEl->nnode != node)
			return 0;
		if(pCacheEl->npars != parsize)
			return 0;
		if(pCacheEl->nPool != poolsize)
			return 0;
		for(i = 0; i < poolsize; i++) {
			if(pCacheEl->pPool[i] != m_parBuff1[i])
				return 0;
		}

		for(i = 0; i < parsize; i++) {
			parset[i] = m_pRorderInverse[pCacheEl->pPars[i] - 1] - 1;
		}
		*probNode = *pCacheEl->pNodeProb;
		*pflik = pCacheEl->fLogLik;
#ifdef DEBUG_INFO
Rprintf("\n    HIT, parset = ");
for(i = 0; i < parsize; i++) {
	sprintf(str,"%d  ", parset[i]);Rprintf(str);
}Rprintf("\n");
#endif
		return 1;
	}

int c_cache::setCachedProb(int *ppool, int poolsize, int node, int *parset, int parsize, 
				PROB_LIST<double> *probNode, double flik) {
		int i, j, prime;
		unsigned nlookup;

		if(!m_bUseCache)
			return 0;
		if(m_numNodes < 1 || !m_pRorder)
			return 0;

		if(!g_pcache) {
			nlookup = m_numNodes;
			for(i = 0; i < m_maxParentSet; i++)
				nlookup *= m_numNodes;
			if(nlookup < PRIMES_1000[PRIMES_NUM-1]) {
				for(i = 0; i < PRIMES_NUM; i++)
					if(nlookup >= PRIMES_1000[i]) {
						nlookup = PRIMES_1000[i];
						break;
					}
			}
			else {
				prime = PRIMES_1000[PRIMES_NUM-1];
				for(i = 0; i < PRIMES_NUM; i++) {
					if(PRIMES_1000[i] >= (unsigned)m_numNodes) {
						prime = PRIMES_1000[i];
						break;
					}
				}
				nlookup = prime;
				for(i = 0; i < m_maxParentSet; i++)
					nlookup *= prime;
			}
			if(nlookup < PRIMES_1000[10])
				nlookup = PRIMES_1000[10];
			
			g_nCacheBits = 1;
			i = 1;
			while(i < m_numNodes*m_maxParentSet) {
				g_nCacheBits++;
				i <<= 1;
			}
#ifdef DEBUG_INFO
sprintf(str,"nCacheBits = %d  ", g_nCacheBits); Rprintf(str);
#endif
			InitializeCache(nlookup, g_nCacheBits);
		}

#ifdef DEBUG_INFO
sprintf(str,"setCachedProb node=%d, poolsize = %d, parsize = %d, pool = ", node, poolsize, parsize); Rprintf(str);
for(i = 0; i < poolsize; i++) {
	sprintf(str,"%d  ", ppool[i]);Rprintf(str);
}Rprintf("\n");
#endif

		node = m_pRorder[node];
		for(i = 0; i < poolsize; i++)
			m_parBuff1[i] = m_pRorder[ppool[i]];
		_quick_sort<int>(m_parBuff1, poolsize);

		for(i = 0; i < parsize; i++)
			m_parBuff2[i] = m_pRorder[parset[i]];

#ifdef DEBUG_INFO
sprintf(str,"    reordered node=%d, pool = ", node);Rprintf(str);
for(i = 0; i < poolsize; i++) {
	sprintf(str,"%d  ", m_parBuff1[i]);Rprintf(str);
}Rprintf("\n par = ");
for(i = 0; i < parsize; i++) {
	sprintf(str,"%d  ", m_parBuff2[i]);Rprintf(str);
}Rprintf("\n");
#endif
		CATNET_CACHE_EL<double> *pCacheEl = new CATNET_CACHE_EL<double>(
			m_parBuff1, poolsize, node, m_parBuff2, parsize, probNode, flik);

		nlookup = 1;
		for(i = 0; i < poolsize; i++) {
			j = m_parBuff1[i] - 1;
			while(j >= PRIMES_NUM)
				j -= PRIMES_NUM;
			nlookup *= PRIMES_1000[PRIMES_NUM - j - 1];
			while(nlookup >= g_ncache)
				nlookup -= g_ncache;		
		}
		nlookup = (nlookup << g_nCacheBits) + node + m_numNodes*parsize;
		while(nlookup >= g_ncache)
			nlookup -= g_ncache;	

#ifdef DEBUG_INFO
sprintf(str,"\nnlookup=%d\n", nlookup);Rprintf(str);
#endif
		if(g_pcache[nlookup]) {

			delete g_pcache[nlookup];
		}
		g_pcache[nlookup] = pCacheEl;
		/*clock_t lt = clock();
		Rprintf("%ld\n", lt);*/
		return 1;
	}

