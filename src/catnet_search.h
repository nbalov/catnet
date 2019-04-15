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
 * catnet_search.h
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
#include "search_params.h"

#ifndef CATNET_SEARCH_H
#define CATNET_SEARCH_H

using namespace std;

template<class t_node, int t_node_size, class t_prob>
class CATNET_SEARCH:  public c_thread, public c_cache {

protected:
	int m_nCatnets;
	CATNET<t_node, t_node_size, t_prob> **m_pCatnets;

	int m_numNodes, *m_pNodeNumCats, **m_pNodeCats, m_numSamples;

public:
	CATNET_SEARCH() {
		m_nCatnets = 0;
		m_pCatnets = 0;
		m_numNodes = 0;
		m_pNodeNumCats = 0;
		m_pNodeCats = 0;
	}

	~CATNET_SEARCH() {
		_release();
	}

	int numCatnets() {
		return m_nCatnets;
	}

	CATNET<t_node, t_node_size, t_prob> **catnets() {
		return m_pCatnets;
	}

protected:
	void _release() {
		int i;
		if(m_pCatnets) {
			for(i = 0; i < m_nCatnets; i++)
				if(m_pCatnets[i]) {
					delete m_pCatnets[i];
					m_pCatnets[i] = 0;
				}
			CATNET_FREE(m_pCatnets);
		}
		m_pCatnets = 0;
		m_nCatnets = 0;
		if(m_pNodeCats) {
			for(i = 0; i < m_numNodes; i++) 
				if(m_pNodeCats[i])
					CATNET_FREE(m_pNodeCats[i]);
			CATNET_FREE(m_pNodeCats);
		}
		if(m_pNodeNumCats) 
			CATNET_FREE(m_pNodeNumCats);
	}

	/* returns increasing subsets of 'parset' of size 'parsize' */
	void combinationSets(int **&plist, int &nlist, int *curset, int *parset, int nparset, int parid, int parsize) {
		int i, ancestor;
 
		if(parid < 0 || parid >= parsize)
			return;
		ancestor = -1;
		if(parid > 0)
			ancestor = curset[parid-1];

		if(parid == parsize - 1) {
			for(i = 0; i < nparset; i++) {
				if(parset[i] <= ancestor)
					continue;
				int **pnewlist = (int**)CATNET_MALLOC((nlist+1)*sizeof(int*));
				if(pnewlist && nlist > 0)
					memcpy(pnewlist, plist, nlist*sizeof(int*));
				pnewlist[nlist] = (int*)CATNET_MALLOC(parsize*sizeof(int));
				if(pnewlist[nlist] && curset) {
					memcpy(pnewlist[nlist], curset, (parsize-1)*sizeof(int));
				}
				pnewlist[nlist][parsize-1] = parset[i];

				CATNET_FREE(plist);
				plist = pnewlist;
				nlist++;
			}
			if(curset) {
				CATNET_FREE(curset);
				curset = 0;
			}
			return;
		}
		for(i = 0; i < nparset; i++) {
			if(parset[i] <= ancestor)
				continue;
			int *pnewset = (int*)CATNET_MALLOC((parid+1)*sizeof(int));
			if(pnewset && curset && parid > 0)
				memcpy(pnewset, curset, parid*sizeof(int));
			pnewset[parid] = parset[i];
			combinationSets(plist, nlist, pnewset, parset, nparset, parid+1, parsize);
		}
		if(curset) {
			CATNET_FREE(curset);
			curset = 0;
		}
	}

public:
	/* psamples and perturbations are sample=columns and node=rows. */
	/* Each parentsPool[i] is numnodes long ! */
	int estimate(SEARCH_PARAMETERS *pestim) {

		if(!pestim)
			return 0;
		int numnodes = pestim->m_numNodes;
		int numsamples = pestim->m_numSamples;
		int *psamples = pestim->m_pSamples;
		int *perturbations = pestim->m_pPerturbations;
		int maxParentSet = pestim->m_maxParentSet;
		int *parSizes = pestim->m_pParentSizes;
		int maxComplexity = pestim->m_maxComplexity;
		int **parentsPool = pestim->m_parentsPool;
		int **fixedParentsPool = pestim->m_fixedParentsPool;
		int becho = pestim->m_echo;

		int i, j, k, d, ncomb, ncombMaxLogLik, nnode, nodecomplx;
		int mincat, maxcat, nocache;
		int maxCategories, numsubsamples, complx, bEqualCategories;
		int *parset, parsetsize, *idparset, *fixparset, fixparsetsize;
		int *paux, *psubsamples, **pcomblist, ncomblist, maxpars, ballow, bfixallow;
		CATNET<t_node, t_node_size, t_prob> baseCatnet, *pNewNet, **pCurCatnetList;

		t_prob fLogLik, fMaxLogLik, tempLogLik;
		PROB_LIST<t_prob> probMaxNode, *pProbNode;

		_release();

		if(numnodes < 1 || numsamples < 1 || !psamples)
			return CATNET_ERR_PARAM;

		if(maxComplexity < numnodes)
			maxComplexity = numnodes;
		
		m_numNodes = numnodes;
		m_numSamples = numsamples;

		maxCategories = 0;

		m_pNodeNumCats = (int*)CATNET_MALLOC(numnodes*sizeof(int));
		if (!m_pNodeNumCats)
			return CATNET_ERR_MEM;
		m_pNodeCats = (int**)CATNET_MALLOC(numnodes*sizeof(int*));
		if (!m_pNodeCats) {
			CATNET_FREE(m_pNodeNumCats);
			return CATNET_ERR_MEM;
		}
		memset(m_pNodeCats,    0, numnodes*sizeof(int*));
		memset(m_pNodeNumCats, 0, numnodes*sizeof(int));

		if(pestim->m_pNodeNumCats && pestim->m_pNodeCats) {
			if (m_pNodeNumCats)
				memcpy(m_pNodeNumCats, pestim->m_pNodeNumCats, numnodes*sizeof(int));
			for(i = 0; i < numnodes; i++) {
				m_pNodeCats[i] = (int*)CATNET_MALLOC(m_pNodeNumCats[i]*sizeof(int));
				if (m_pNodeCats[i] && pestim->m_pNodeCats[i])
					memcpy(m_pNodeCats[i], pestim->m_pNodeCats[i], m_pNodeNumCats[i]*sizeof(int));
			}
		}
		else { 
			for(i = 0; i < numnodes; i++) {
				mincat = INT_MAX;
				maxcat = -INT_MAX;
				for(j = 0; j < numsamples; j++) {
					if(psamples[j*numnodes + i] == CATNET_NAN)
						continue;
					if(psamples[j*numnodes + i] < mincat)
						mincat = psamples[j*numnodes + i];
					if(psamples[j*numnodes + i] > maxcat)
						maxcat = psamples[j*numnodes + i];
				}
				m_pNodeNumCats[i] = maxcat - mincat + 1;
				m_pNodeCats[i] = (int*)CATNET_MALLOC(m_pNodeNumCats[i]*sizeof(int));
				for(j = 0; j < m_pNodeNumCats[i]; j++)
					m_pNodeCats[i][j] = mincat + j;
				/* order m_pNodeNumCats[i] */
				for(j = 0; j < m_pNodeNumCats[i]; j++) {
					for(k = j + 1; k < m_pNodeNumCats[i]; k++) {
						if(m_pNodeCats[i][j] > m_pNodeCats[i][k]) {
							d = m_pNodeCats[i][j]; 
							m_pNodeCats[i][j] = m_pNodeCats[i][k];
							m_pNodeCats[i][k] = d;
						}
					}
				}
			}
		}
		for(i = 0; i < numnodes; i++) {
			for(j = 0; j < numsamples; j++) {
				if(psamples[j*numnodes + i] == CATNET_NAN) {
					// CATNET_NANs will be assigned m_pNodeNumCats[i]
					psamples[j*numnodes + i] = m_pNodeNumCats[i];
					continue;
				}
				for(d = 0; d < m_pNodeNumCats[i]; d++)
					if(m_pNodeCats[i][d] == psamples[j*numnodes + i])
						break;
				if(d >= m_pNodeNumCats[i])
					return CATNET_ERR_PARAM;
				psamples[j*numnodes + i] = d;
			}
			if(maxCategories < m_pNodeNumCats[i])
				maxCategories = m_pNodeNumCats[i];
			if(i > 1 && m_pNodeNumCats[i] != m_pNodeNumCats[0])
				bEqualCategories = 0;
		}


		bEqualCategories = 1;
		for(i = 0; i < numnodes; i++) {
			if(i > 1 && m_pNodeNumCats[i] != m_pNodeNumCats[0])
				bEqualCategories = 0;
		}

		parset    = (int*)CATNET_MALLOC(numnodes*sizeof(int));
		idparset  = (int*)CATNET_MALLOC(numnodes*sizeof(int));
		fixparset = (int*)CATNET_MALLOC(numnodes*sizeof(int));

		m_nCatnets = maxComplexity + 1;
		m_pCatnets = (CATNET<t_node, t_node_size, t_prob>**)CATNET_MALLOC(m_nCatnets*sizeof(CATNET<t_node, t_node_size, t_prob>*));

		pCurCatnetList = (CATNET<t_node, t_node_size, t_prob>**)CATNET_MALLOC(m_nCatnets*sizeof(CATNET<t_node, t_node_size, t_prob>*));

		psubsamples = 0;
		if(perturbations) {
			psubsamples = (int*)CATNET_MALLOC(numnodes*numsamples*sizeof(int));
		}

		/* create a network without edges*/
		pNewNet = new CATNET<t_node, t_node_size, t_prob>
				(numnodes, 0/*maxParentSet*/, maxCategories, 0, 0, 0, m_pNodeNumCats);

		if (!parset || !idparset || !fixparset || !m_pCatnets || !pNewNet || !pCurCatnetList) {
			if (parset) 
				CATNET_FREE(parset);
			if (idparset) 
				CATNET_FREE(idparset);
			if (fixparset) 
				CATNET_FREE(fixparset);
			if (m_pNodeCats)
				CATNET_FREE(m_pNodeCats);
			m_pNodeCats = 0;
			if (m_pNodeNumCats)
				CATNET_FREE(m_pNodeNumCats);
			m_pNodeNumCats = 0;
			if (m_pCatnets)
				CATNET_FREE(m_pCatnets);
			m_pCatnets = 0;
			if (pCurCatnetList)
				CATNET_FREE(pCurCatnetList);
			if (psubsamples)
				CATNET_FREE(psubsamples);
			if (pNewNet)
				delete pNewNet;
			return CATNET_ERR_MEM;
		}

		memset(m_pCatnets, 0, m_nCatnets*sizeof(CATNET<t_node, t_node_size, t_prob>*));

		/* set parents */
		for(nnode = 0; nnode < numnodes; nnode++) {
			if(fixedParentsPool && fixedParentsPool[nnode]) {
				fixparsetsize = 0;
				for(j = 0; j < nnode; j++) {
					ballow = 1;
					if(parentsPool && parentsPool[nnode]) {
						ballow = 0;
						for(k = 0; k < numnodes; k++) {
							if(j == parentsPool[nnode][k])
								ballow = 1;
						}
					}
					if(parentsPool && !parentsPool[nnode])
						ballow = 0;
					bfixallow = 0;
					if(fixedParentsPool && fixedParentsPool[nnode]) {
						for(k = 0; k < numnodes; k++) {
							if(j == fixedParentsPool[nnode][k])
								bfixallow = 1;
						}
					}
					if(!ballow)
						continue;
					if(bfixallow) {
					  fixparset[fixparsetsize] = j;
					  fixparsetsize++;
					}
				}
				if(fixparsetsize > 0)
					pNewNet -> setParents(nnode, fixparset, fixparsetsize);
			}
		}

		baseCatnet.init(numnodes, maxParentSet, maxCategories, 0, 0, 0, m_pNodeNumCats);

		// set sample probabilities and calculate log-likelihood
		for(nnode = 0; nnode < numnodes; nnode++) {
			if(perturbations && psubsamples) {
				numsubsamples = 0;
				for(j = 0; j < numsamples; j++) {
					if(!perturbations[j * numnodes + nnode]) {
						memcpy(psubsamples + numsubsamples*numnodes, psamples + j*numnodes, numnodes*sizeof(int));
						numsubsamples++;
					}
				}
				pNewNet->setNodeSampleProb(nnode, psubsamples, numsubsamples);
			}
			else {
				pNewNet->setNodeSampleProb(nnode, psamples, numsamples);
			}
		}

		complx = pNewNet->complexity();
		m_pCatnets[complx] = pNewNet;

		/* main loop of consequential non-empty-parenthood-node additions */
		for(nnode = 1; nnode < numnodes; nnode++) {

			if(_wait_stop_event(4/*millisecs*/) == 0)
				break;

			if(becho) {
				Rprintf("processing node %d\n", nnode+1);
				Rprintf("    [#parents][#combinations] = ");
			}

			fixparsetsize = 0;
			parsetsize = 0;
			for(j = 0; j < nnode; j++) {
				ballow = 1;
				if(parentsPool && parentsPool[nnode]) {
					ballow = 0;
					for(k = 0; k < numnodes; k++) {
						if(j == parentsPool[nnode][k])
							ballow = 1;
					}
				}
				if(parentsPool && !parentsPool[nnode])
					ballow = 0;
				bfixallow = 0;
				if(fixedParentsPool && fixedParentsPool[nnode]) {
					for(k = 0; k < numnodes; k++) {
						if(j == fixedParentsPool[nnode][k])
							bfixallow = 1;
					}
				}
				if(!ballow)
					continue;
				if(bfixallow) {
				  fixparset[fixparsetsize] = j;
				  fixparsetsize++;
				}
				else {
				  parset[parsetsize] = j;
				  parsetsize++;
				}
			}
			/* extend the content before sending to cache; parsetsize + fixparsetsize < numnodes */
			if (parset && fixparset && fixparsetsize > 0)
				memcpy(parset + parsetsize, fixparset, fixparsetsize*sizeof(int));

			/* check out wheather the parent pool has equal number of categories */
			bEqualCategories = 1;
			for(j = 0; j < parsetsize + fixparsetsize; j++) 
				if(j > 0 && m_pNodeNumCats[parset[j]] != m_pNodeNumCats[parset[0]])
					bEqualCategories = 0;

			maxpars = maxParentSet;
			if(parSizes && parSizes[nnode] < maxParentSet)
				maxpars = parSizes[nnode];

			if(maxpars > parsetsize + fixparsetsize)
				maxpars = parsetsize + fixparsetsize;

			memset(pCurCatnetList, 0, m_nCatnets*sizeof(CATNET<t_node, t_node_size, t_prob>*));

			if(bEqualCategories) {

			for(d = fixparsetsize + 1; d <= maxpars; d++) {

				//if(_wait_stop_event(1/*millisecs*/) == 0)
				//	break;

				nocache = 1;
				if(!pestim->m_pCacheMutex) {
					nocache = !getCachedProb(parset, parsetsize + fixparsetsize, nnode, 
						idparset, d, &probMaxNode, &fMaxLogLik);
				}
				else {
					MUTEX_LOCK(pestim->m_pCacheMutex);
					nocache = !getCachedProb(parset, parsetsize + fixparsetsize, nnode, 
						idparset, d, &probMaxNode, &fMaxLogLik);
					MUTEX_UNLOCK(pestim->m_pCacheMutex);
				}

				if(nocache) { 

					pcomblist = 0;
					ncomblist = 0;
					combinationSets(pcomblist, ncomblist, 0, parset, parsetsize, 0, d - fixparsetsize);

				        if(fixparsetsize > 0) {
				        	if(!pcomblist || ncomblist < 1) {
				        	    	pcomblist = (int**)CATNET_MALLOC(1*sizeof(int*));
							if (!pcomblist)
								return CATNET_ERR_MEM;
				            		pcomblist[0] = 0;	
				            		ncomblist = 1;
				          	}
				        	for(k = 0; k < ncomblist; k++) {
				            		paux = (int*)CATNET_MALLOC(d*sizeof(int));
							if (!paux)
								return CATNET_ERR_MEM;
							for(j = 0; j < fixparsetsize; j++)
				            			paux[j] = fixparset[j];
					        	if(paux && pcomblist[k] && d > fixparsetsize) {
				            			memcpy(paux + fixparsetsize, pcomblist[k], (d-fixparsetsize)*sizeof(int));
							}
				            		if(pcomblist[k])
				            			CATNET_FREE(pcomblist[k]); 
				           		pcomblist[k] = paux;
						}
					}
			
					if(becho)
						Rprintf("[%d]%d  ", d, ncomblist);

					fMaxLogLik = -FLT_MAX;
					ncombMaxLogLik = -1;
					probMaxNode.reset();

					for(ncomb = 0; ncomb < ncomblist; ncomb++) {
						// add pcomplist[j] parent set to nnode
						baseCatnet.setParents(nnode, pcomblist[ncomb], d);
			     
						// add perturbation
						if(perturbations && psubsamples) {
							numsubsamples = 0;
							for(j = 0; j < numsamples; j++) {
								if(!perturbations[j * numnodes + nnode]) {
									memcpy(psubsamples + numsubsamples*numnodes, psamples + j*numnodes, numnodes*sizeof(int));
									numsubsamples++;
								}
							}
							fLogLik = baseCatnet.setNodeSampleProb(nnode, psubsamples, numsubsamples);
						}
						else {
							fLogLik = baseCatnet.setNodeSampleProb(nnode, psamples, numsamples);
						}

						if(fMaxLogLik < fLogLik) {
							pProbNode = baseCatnet.getNodeProb(nnode); 
							fMaxLogLik = fLogLik;
							ncombMaxLogLik = ncomb;
							if(pProbNode)
								probMaxNode = *pProbNode;
						}
					} /* for ncomb */

					if(idparset && ncombMaxLogLik >= 0 && d > 0)
						memcpy(idparset, pcomblist[ncombMaxLogLik], d*sizeof(int));

					/* release combination set */
        				for(ncomb = 0; ncomb < ncomblist; ncomb++) {
        	  				if(pcomblist[ncomb])
        	    					CATNET_FREE(pcomblist[ncomb]);
        	  				pcomblist[ncomb] = NULL;
					} /* for ncomb */
        				CATNET_FREE(pcomblist);
        				pcomblist = 0;
					ncomblist = 0;

					if(ncombMaxLogLik < 0){
						/* retain the same m_pCatnets list */
						continue;
					}

					if(pestim->m_pCacheMutex) {
						MUTEX_LOCK(pestim->m_pCacheMutex);
						setCachedProb(parset, parsetsize + fixparsetsize, 
							nnode, idparset, d, &probMaxNode, fMaxLogLik);
						MUTEX_UNLOCK(pestim->m_pCacheMutex);
					}
					else
						setCachedProb(parset, parsetsize + fixparsetsize, 
							nnode, idparset, d, &probMaxNode, fMaxLogLik);
	
				} /* if(!getCachedProb) */

				nodecomplx = m_pNodeNumCats[nnode]-1;
				for(k = 0; k < d; k++)
					nodecomplx *= m_pNodeNumCats[idparset[k]];

				for(k = 0; k < m_nCatnets; k++) {
					if(!m_pCatnets[k])
						continue;
					complx = m_pCatnets[k]->complexity() - 
						m_pCatnets[k]->nodeComplexity(nnode) +
						nodecomplx;
					pProbNode = m_pCatnets[k]->getNodeProb(nnode);
					if(complx > maxComplexity || !pProbNode) {
						continue;
					}
					tempLogLik = m_pCatnets[k]->loglik() - pProbNode->loglik + fMaxLogLik;
					if(!pCurCatnetList[complx] && tempLogLik > -FLT_MAX) {
						pCurCatnetList[complx] = new CATNET<t_node, t_node_size, t_prob>;
					}
					if(pCurCatnetList[complx] && 
						pCurCatnetList[complx]->loglik() < tempLogLik) {
							*pCurCatnetList[complx] = *m_pCatnets[k];
							pCurCatnetList[complx]->setParents(nnode, idparset, d);
							pCurCatnetList[complx]->setNodeProb(nnode, &probMaxNode);
						}
					}
				} /* for k */

			} /* for d */

			else /*if(!bEqualCategories)*/ {

			for(d = fixparsetsize + 1; d <= maxpars; d++) {
				
				if(_wait_stop_event(0/*millisecs*/) == 0)
					break;

				pcomblist = 0;
				ncomblist = 0;
				combinationSets(pcomblist, ncomblist, 0, parset, parsetsize, 0, d - fixparsetsize);

				if(fixparsetsize > 0) {
					if(!pcomblist || ncomblist < 1) {
					    	pcomblist = (int**)CATNET_MALLOC(1*sizeof(int*));
						if (!pcomblist)
							return CATNET_ERR_MEM;
			         		pcomblist[0] = 0;	
						ncomblist = 1;
					}
					for(k = 0; k < ncomblist; k++) {
				        	paux = (int*)CATNET_MALLOC(d*sizeof(int));
						if (!paux)
							return CATNET_ERR_MEM;
						for(j = 0; j < fixparsetsize; j++)
				            		paux[j] = fixparset[j];
					        if(paux && pcomblist[k] && d > fixparsetsize) {
				            		memcpy(paux + fixparsetsize, pcomblist[k], (d-fixparsetsize)*sizeof(int));
						}
				            	if(pcomblist[k])
				            		CATNET_FREE(pcomblist[k]); 
				           	pcomblist[k] = paux;
					}
				}
			
				if(becho)
					Rprintf("[%d]%d  ", d, ncomblist);

				fMaxLogLik = -FLT_MAX;
				ncombMaxLogLik = -1;
				probMaxNode.reset();

				for(ncomb = 0; ncomb < ncomblist; ncomb++) {

					nocache = 1;
					if(!pestim->m_pCacheMutex) {
						nocache = !getCachedProb(pcomblist[ncomb], d, nnode, idparset, d, &probMaxNode, &fMaxLogLik);
					}
					else {
						MUTEX_LOCK(pestim->m_pCacheMutex);
						nocache = !getCachedProb(pcomblist[ncomb], d, nnode, idparset, d, &probMaxNode, &fMaxLogLik);
						MUTEX_UNLOCK(pestim->m_pCacheMutex);
					}

					if(nocache) { 
						if (idparset && pcomblist[ncomb] && d > 0)
							memcpy(idparset, pcomblist[ncomb], d*sizeof(int));
						// add pcomplist[j] parent set to nnode
						baseCatnet.setParents(nnode, idparset, d);
			     
						// add perturbation
						if(perturbations && psubsamples) {
							numsubsamples = 0;
							for(j = 0; j < numsamples; j++) {
								if(!perturbations[j * numnodes + nnode]) {
									memcpy(psubsamples + numsubsamples*numnodes, psamples + j*numnodes, numnodes*sizeof(int));
									numsubsamples++;
								}
							}
							fMaxLogLik = baseCatnet.setNodeSampleProb(nnode, psubsamples, numsubsamples);
						}
						else {
							fMaxLogLik = baseCatnet.setNodeSampleProb(nnode, psamples, numsamples);
						}

						pProbNode = baseCatnet.getNodeProb(nnode);
						if(pProbNode)
							probMaxNode = *pProbNode;
						
						if(pestim->m_pCacheMutex) {
							MUTEX_LOCK(pestim->m_pCacheMutex);
							setCachedProb(idparset, d, nnode, idparset, d,
									&probMaxNode, fMaxLogLik);
							MUTEX_UNLOCK(pestim->m_pCacheMutex);
						}
						else
							setCachedProb(idparset, d, nnode, idparset, d,
									&probMaxNode, fMaxLogLik);

					} /* if(!getCachedProb) */

					/* find nnode-complexity for pcomblist[ncomb] parent set*/
					nodecomplx = m_pNodeNumCats[nnode]-1;
					for(k = 0; k < d; k++)
						nodecomplx *= m_pNodeNumCats[idparset[k]];
					for(k = 0; k < m_nCatnets; k++) {
						if(!m_pCatnets[k])
							continue;
						pProbNode = m_pCatnets[k]->getNodeProb(nnode);
						complx = m_pCatnets[k]->complexity() - 
							m_pCatnets[k]->nodeComplexity(nnode) +
							nodecomplx;
						if(complx > maxComplexity || !pProbNode) 
							continue;
						tempLogLik = m_pCatnets[k]->loglik() - 
							pProbNode->loglik + fMaxLogLik;
						if(!pCurCatnetList[complx] && tempLogLik > -FLT_MAX) {
							pCurCatnetList[complx] = new CATNET<t_node, t_node_size, t_prob>;
							if (!pCurCatnetList[complx])
								return CATNET_ERR_MEM;
						}
						if(pCurCatnetList[complx] && 
							pCurCatnetList[complx]->loglik() < tempLogLik) {
								*pCurCatnetList[complx] = *m_pCatnets[k];
								pCurCatnetList[complx]->setParents(nnode, idparset, d);
								pCurCatnetList[complx]->setNodeProb(nnode, &probMaxNode);
						}
					} /* for k */

				} /* for ncomb */
	
				/* release combination set */
        			for(ncomb = 0; ncomb < ncomblist; ncomb++) {
          				if(pcomblist[ncomb])
            					CATNET_FREE(pcomblist[ncomb]);
          				pcomblist[ncomb] = NULL;
				} /* for ncomb */
        			CATNET_FREE(pcomblist);
        			pcomblist = 0;
				ncomblist = 0;

			} /* for d */

			} /* if(!bEqualCategories) */

			if(becho)
				Rprintf("\n");

			for(j = 0; j < m_nCatnets; j++) {
				if(m_pCatnets[j]) {
					if(pCurCatnetList[j]) {
						if(m_pCatnets[j]->loglik() < pCurCatnetList[j]->loglik()) {
							delete m_pCatnets[j];
							m_pCatnets[j] = pCurCatnetList[j];
							pCurCatnetList[j] = 0;
						}
						else {
							delete pCurCatnetList[j];
							pCurCatnetList[j] = 0;
						}
					}
				}
				else {
					m_pCatnets[j] = pCurCatnetList[j];
					pCurCatnetList[j] = 0;
				}
			}
		} // for(nnode = 0; nnode < numnodes; nnode++)

		CATNET_FREE(pCurCatnetList);
		CATNET_FREE(parset);
		CATNET_FREE(fixparset);
		CATNET_FREE(idparset);

		if(psubsamples)
			CATNET_FREE(psubsamples);

		for(j = 0; j < m_nCatnets; j++) {
			if(m_pCatnets[j]) {
				m_pCatnets[j] -> normalizeProbabilities();
				m_pCatnets[j] -> setCategoryIndices(m_pNodeNumCats, m_pNodeCats);
			}
		} 

		if(m_pNodeCats) {
			for(i = 0; i < m_numNodes; i++) 
				if(m_pNodeCats[i])
					CATNET_FREE(m_pNodeCats[i]);
			CATNET_FREE(m_pNodeCats);
		}
		m_pNodeCats = 0;
		if(m_pNodeNumCats) 
			CATNET_FREE(m_pNodeNumCats);
		m_pNodeNumCats = 0;
		
		return 0;
	}

	/* psamples and perturbations are sample=columns and node=rows. */
	/* Each parentsPool[i] is numnodes long ! */
	int pairwise_KL_distance(SEARCH_PARAMETERS *pestim) {

		if(!pestim)
			return 0;
		int numnodes = pestim->m_numNodes;
		int numsamples = pestim->m_numSamples;
		int *psamples = pestim->m_pSamples;
		int *perturbations = pestim->m_pPerturbations;
		int maxParentSet = pestim->m_maxParentSet;
		int *parSizes = pestim->m_pParentSizes;
		int maxComplexity = pestim->m_maxComplexity;
		int **parentsPool = pestim->m_parentsPool;
		int **fixedParentsPool = pestim->m_fixedParentsPool;
		int becho = pestim->m_echo;

		int i, j, k, d, ncomb, ncombMaxLogLik, nnode, nodecomplx;
		int mincat, maxcat, nocache;
		int maxCategories, numsubsamples, complx, bEqualCategories;
		int *parset, parsetsize, *idparset, *fixparset, fixparsetsize;
		int *paux, *psubsamples, **pcomblist, ncomblist, maxpars, ballow, bfixallow;
		CATNET<t_node, t_node_size, t_prob> baseCatnet, *pNewNet, **pCurCatnetList;

		t_prob fLogLik, fMaxLogLik, tempLogLik;
		PROB_LIST<t_prob> probMaxNode, *pProbNode;

		_release();

		if(numnodes < 1 || numsamples < 1 || !psamples)
			return 0;
		if(maxComplexity < numnodes)
			maxComplexity = numnodes;
		
		m_numNodes = numnodes;
		m_numSamples = numsamples;

		maxCategories = 0;

		m_pNodeNumCats = (int*)CATNET_MALLOC(numnodes*sizeof(int));
		if (!m_pNodeNumCats)
			return CATNET_ERR_MEM;
		m_pNodeCats = (int**)CATNET_MALLOC(numnodes*sizeof(int*));
		if (!m_pNodeCats) {
			CATNET_FREE(m_pNodeNumCats);
			return CATNET_ERR_MEM;
		}
		memset(m_pNodeCats,    0, numnodes*sizeof(int*));
		memset(m_pNodeNumCats, 0, numnodes*sizeof(int));

		for(i = 0; i < numnodes; i++) {
			mincat = INT_MAX;
			maxcat = -INT_MAX;
			for(j = 0; j < numsamples; j++) {
				if(psamples[j*numnodes + i] < mincat)
					mincat = psamples[j*numnodes + i];
				if(psamples[j*numnodes + i] > maxcat)
					maxcat = psamples[j*numnodes + i];
			}
			m_pNodeNumCats[i] = maxcat - mincat + 1;
			m_pNodeCats[i] = (int*)CATNET_MALLOC(m_pNodeNumCats[i]*sizeof(int));
			for(j = 0; j < m_pNodeNumCats[i]; j++)
				m_pNodeCats[i][j] = mincat + j;
		}
		for(i = 0; i < numnodes; i++) {
			/* order m_pNodeNumCats[i] */
			for(j = 0; j < m_pNodeNumCats[i]; j++) {
				for(k = j + 1; k < m_pNodeNumCats[i]; k++) {
					if(m_pNodeCats[i][j] > m_pNodeCats[i][k]) {
						d = m_pNodeCats[i][j]; 
						m_pNodeCats[i][j] = m_pNodeCats[i][k];
						m_pNodeCats[i][k] = d;
					}
				}
			} 
			for(j = 0; j < numsamples; j++) {
				for(d = 0; d < m_pNodeNumCats[i]; d++)
					if(m_pNodeCats[i][d] == psamples[j*numnodes + i])
						break;
				psamples[j*numnodes + i] = d;
			}
			if(maxCategories < m_pNodeNumCats[i])
				maxCategories = m_pNodeNumCats[i];
			if(i > 1 && m_pNodeNumCats[i] != m_pNodeNumCats[0])
				bEqualCategories = 0;
		}


		bEqualCategories = 1;
		for(i = 0; i < numnodes; i++) {
			if(i > 1 && m_pNodeNumCats[i] != m_pNodeNumCats[0])
				bEqualCategories = 0;
		}

		parset = (int*)CATNET_MALLOC(numnodes*sizeof(int));
		idparset = (int*)CATNET_MALLOC(numnodes*sizeof(int));
		fixparset = (int*)CATNET_MALLOC(numnodes*sizeof(int));

		m_nCatnets = maxComplexity + 1;
		m_pCatnets = (CATNET<t_node, t_node_size, t_prob>**)CATNET_MALLOC(m_nCatnets*sizeof(CATNET<t_node, t_node_size, t_prob>*));

		pCurCatnetList = (CATNET<t_node, t_node_size, t_prob>**)CATNET_MALLOC(m_nCatnets*sizeof(CATNET<t_node, t_node_size, t_prob>*));

		psubsamples = 0;
		if(perturbations) {
			psubsamples = (int*)CATNET_MALLOC(numnodes*numsamples*sizeof(int));
		}

		/* create a network without edges*/
		pNewNet = new CATNET<t_node, t_node_size, t_prob>
				(numnodes, 0/*maxParentSet*/, maxCategories, 0, 0, 0, m_pNodeNumCats);

		if (!parset || !idparset || !fixparset || !m_pCatnets || !pNewNet || !pCurCatnetList) {
			if (parset) 
				CATNET_FREE(parset);
			if (idparset) 
				CATNET_FREE(idparset);
			if (fixparset) 
				CATNET_FREE(fixparset);
			if (m_pNodeCats)
				CATNET_FREE(m_pNodeCats);
			m_pNodeCats = 0;
			if (m_pNodeNumCats)
				CATNET_FREE(m_pNodeNumCats);
			m_pNodeNumCats = 0;
			if (m_pCatnets)
				CATNET_FREE(m_pCatnets);
			m_pCatnets = 0;
			if (pCurCatnetList)
				CATNET_FREE(pCurCatnetList);
			if (psubsamples)
				CATNET_FREE(psubsamples);
			if (pNewNet)
				delete pNewNet;
			return CATNET_ERR_MEM;
		}

		memset(m_pCatnets, 0, m_nCatnets*sizeof(CATNET<t_node, t_node_size, t_prob>*));

		/* set parents */
		for(nnode = 0; nnode < numnodes; nnode++) {
			if(fixedParentsPool && fixedParentsPool[nnode]) {
				fixparsetsize = 0;
				for(j = 0; j < nnode; j++) {
					ballow = 1;
					if(parentsPool && parentsPool[nnode]) {
						ballow = 0;
						for(k = 0; k < numnodes; k++) {
							if(j == parentsPool[nnode][k])
								ballow = 1;
						}
					}
					if(parentsPool && !parentsPool[nnode])
						ballow = 0;
					bfixallow = 0;
					if(fixedParentsPool && fixedParentsPool[nnode]) {
						for(k = 0; k < numnodes; k++) {
							if(j == fixedParentsPool[nnode][k])
								bfixallow = 1;
						}
					}
					if(!ballow)
						continue;
					if(bfixallow) {
					  fixparset[fixparsetsize] = j;
					  fixparsetsize++;
					}
				}
				if(fixparsetsize > 0)
					pNewNet -> setParents(nnode, fixparset, fixparsetsize);
			}
		}

		baseCatnet.init(numnodes, maxParentSet, maxCategories, 0, 0, 0, m_pNodeNumCats);

		// set sample probabilities and calculate log-likelihood
		for(nnode = 0; nnode < numnodes; nnode++) {
			if(perturbations && psubsamples) {
				numsubsamples = 0;
				for(j = 0; j < numsamples; j++) {
					if(!perturbations[j * numnodes + nnode]) {
						memcpy(psubsamples + numsubsamples*numnodes, psamples + j*numnodes, numnodes*sizeof(int));
						numsubsamples++;
					}
				}
				pNewNet->setNodeSampleProb(nnode, psubsamples, numsubsamples);
			}
			else {
				pNewNet->setNodeSampleProb(nnode, psamples, numsamples);
			}
		}

		complx = pNewNet->complexity();
		m_pCatnets[complx] = pNewNet;

		/* main loop of consequential non-empty-parenthood-node additions */
		for(nnode = 1; nnode < numnodes; nnode++) {

			if(_wait_stop_event(4/*millisecs*/) == 0)
				break;

			if(becho) {
				Rprintf("processing node %d\n", nnode+1);
				Rprintf("    [#parents][#combinations] = ");
			}

			fixparsetsize = 0;
			parsetsize = 0;
			for(j = 0; j < nnode; j++) {
				ballow = 1;
				if(parentsPool && parentsPool[nnode]) {
					ballow = 0;
					for(k = 0; k < numnodes; k++) {
						if(j == parentsPool[nnode][k])
							ballow = 1;
					}
				}
				if(parentsPool && !parentsPool[nnode])
					ballow = 0;
				bfixallow = 0;
				if(fixedParentsPool && fixedParentsPool[nnode]) {
					for(k = 0; k < numnodes; k++) {
						if(j == fixedParentsPool[nnode][k])
							bfixallow = 1;
					}
				}
				if(!ballow)
					continue;
				if(bfixallow) {
				  fixparset[fixparsetsize] = j;
				  fixparsetsize++;
				}
				else {
				  parset[parsetsize] = j;
				  parsetsize++;
				}
			}
			/* extend the content before sending to cache; parsetsize + fixparsetsize < numnodes */
			if (parset && fixparset && fixparsetsize > 0)
				memcpy(parset + parsetsize, fixparset, fixparsetsize*sizeof(int));

			/* check out wheather the parent pool has equal number of categories */
			bEqualCategories = 1;
			for(j = 0; j < parsetsize + fixparsetsize; j++) 
				if(j > 0 && m_pNodeNumCats[parset[j]] != m_pNodeNumCats[parset[0]])
					bEqualCategories = 0;

			maxpars = maxParentSet;
			if(parSizes && parSizes[nnode] < maxParentSet)
				maxpars = parSizes[nnode];

			if(maxpars > parsetsize + fixparsetsize)
				maxpars = parsetsize + fixparsetsize;

			memset(pCurCatnetList, 0, m_nCatnets*sizeof(CATNET<t_node, t_node_size, t_prob>*));

			if(bEqualCategories) {

			for(d = fixparsetsize + 1; d <= maxpars; d++) {

				//if(_wait_stop_event(1/*millisecs*/) == 0)
				//	break;

				nocache = 1;
				if(!pestim->m_pCacheMutex) {
					nocache = !getCachedProb(parset, parsetsize + fixparsetsize, nnode, 
						idparset, d, &probMaxNode, &fMaxLogLik);
				}
				else {
					MUTEX_LOCK(pestim->m_pCacheMutex);
					nocache = !getCachedProb(parset, parsetsize + fixparsetsize, nnode, 
						idparset, d, &probMaxNode, &fMaxLogLik);
					MUTEX_UNLOCK(pestim->m_pCacheMutex);
				}

				if(nocache) { 

					pcomblist = 0;
					ncomblist = 0;
					combinationSets(pcomblist, ncomblist, 0, parset, parsetsize, 0, d - fixparsetsize);

				        if(fixparsetsize > 0) {
				        	if(!pcomblist || ncomblist < 1) {
				        	    	pcomblist = (int**)CATNET_MALLOC(1*sizeof(int*));
							if (!pcomblist)
								return CATNET_ERR_MEM;
				            		pcomblist[0] = 0;	
				            		ncomblist = 1;
				          	}
				        	for(k = 0; k < ncomblist; k++) {
				            		paux = (int*)CATNET_MALLOC(d*sizeof(int));
							if (!paux)
								return CATNET_ERR_MEM;
							for(j = 0; j < fixparsetsize; j++)
				            			paux[j] = fixparset[j];
					        	if(paux && pcomblist[k] && d > fixparsetsize) {
				            			memcpy(paux + fixparsetsize, pcomblist[k], (d-fixparsetsize)*sizeof(int));
							}
				            		if(pcomblist[k])
				            			CATNET_FREE(pcomblist[k]); 
				           		pcomblist[k] = paux;
						}
					}
			
					if(becho)
						Rprintf("[%d]%d  ", d, ncomblist);

					fMaxLogLik = -FLT_MAX;
					ncombMaxLogLik = -1;
					probMaxNode.reset();

					for(ncomb = 0; ncomb < ncomblist; ncomb++) {
						// add pcomplist[j] parent set to nnode
						baseCatnet.setParents(nnode, pcomblist[ncomb], d);
			     
						// add perturbation
						if(perturbations && psubsamples) {
							numsubsamples = 0;
							for(j = 0; j < numsamples; j++) {
								if(!perturbations[j * numnodes + nnode]) {
									memcpy(psubsamples + numsubsamples*numnodes, psamples + j*numnodes, numnodes*sizeof(int));
									numsubsamples++;
								}
							}
							fLogLik = baseCatnet.setNodeSampleProb(nnode, psubsamples, numsubsamples);
						}
						else {
							fLogLik = baseCatnet.setNodeSampleProb(nnode, psamples, numsamples);
						}

						if(fMaxLogLik < fLogLik) {
							pProbNode = baseCatnet.getNodeProb(nnode); 
							fMaxLogLik = fLogLik;
							ncombMaxLogLik = ncomb;
							if(pProbNode)
								probMaxNode = *pProbNode;
						}
					} /* for ncomb */

					if(idparset && ncombMaxLogLik >= 0 && d > 0)
						memcpy(idparset, pcomblist[ncombMaxLogLik], d*sizeof(int));

					/* release combination set */
        				for(ncomb = 0; ncomb < ncomblist; ncomb++) {
        	  				if(pcomblist[ncomb])
        	    					CATNET_FREE(pcomblist[ncomb]);
        	  				pcomblist[ncomb] = NULL;
					} /* for ncomb */
        				CATNET_FREE(pcomblist);
        				pcomblist = 0;
					ncomblist = 0;

					if(ncombMaxLogLik < 0){
						/* retain the same m_pCatnets list */
						continue;
					}

					if(pestim->m_pCacheMutex) {
						MUTEX_LOCK(pestim->m_pCacheMutex);
						setCachedProb(parset, parsetsize + fixparsetsize, 
							nnode, idparset, d, &probMaxNode, fMaxLogLik);
						MUTEX_UNLOCK(pestim->m_pCacheMutex);
					}
					else
						setCachedProb(parset, parsetsize + fixparsetsize, 
							nnode, idparset, d, &probMaxNode, fMaxLogLik);
	
				} /* if(!getCachedProb) */

				nodecomplx = m_pNodeNumCats[nnode]-1;
				for(k = 0; k < d; k++)
					nodecomplx *= m_pNodeNumCats[idparset[k]];

				for(k = 0; k < m_nCatnets; k++) {
					if(!m_pCatnets[k])
						continue;
					complx = m_pCatnets[k]->complexity() - 
						m_pCatnets[k]->nodeComplexity(nnode) +
						nodecomplx;
					pProbNode = m_pCatnets[k]->getNodeProb(nnode);
					if(complx > maxComplexity || !pProbNode) {
						continue;
					}
					tempLogLik = m_pCatnets[k]->loglik() - pProbNode->loglik + fMaxLogLik;
					if(!pCurCatnetList[complx] && tempLogLik > -FLT_MAX) {
						pCurCatnetList[complx] = new CATNET<t_node, t_node_size, t_prob>;
					}
					if(pCurCatnetList[complx] && 
						pCurCatnetList[complx]->loglik() < tempLogLik) {
							*pCurCatnetList[complx] = *m_pCatnets[k];
							pCurCatnetList[complx]->setParents(nnode, idparset, d);
							pCurCatnetList[complx]->setNodeProb(nnode, &probMaxNode);
						}
					}
				} /* for k */

			} /* for d */

			else /*if(!bEqualCategories)*/ {

			for(d = fixparsetsize + 1; d <= maxpars; d++) {
				
				if(_wait_stop_event(0/*millisecs*/) == 0)
					break;

				pcomblist = 0;
				ncomblist = 0;
				combinationSets(pcomblist, ncomblist, 0, parset, parsetsize, 0, d - fixparsetsize);

				if(fixparsetsize > 0) {
					if(!pcomblist || ncomblist < 1) {
					    	pcomblist = (int**)CATNET_MALLOC(1*sizeof(int*));
						if (!pcomblist)
							return CATNET_ERR_MEM;
			         		pcomblist[0] = 0;	
						ncomblist = 1;
					}
					for(k = 0; k < ncomblist; k++) {
				        	paux = (int*)CATNET_MALLOC(d*sizeof(int));
						if (!paux)
							return CATNET_ERR_MEM;
						for(j = 0; j < fixparsetsize; j++)
				            		paux[j] = fixparset[j];
					        if(paux && pcomblist[k] && d > fixparsetsize) {
				            		memcpy(paux + fixparsetsize, pcomblist[k], (d-fixparsetsize)*sizeof(int));
						}
				            	if(pcomblist[k])
				            		CATNET_FREE(pcomblist[k]); 
				           	pcomblist[k] = paux;
					}
				}
			
				if(becho)
					Rprintf("[%d]%d  ", d, ncomblist);

				fMaxLogLik = -FLT_MAX;
				ncombMaxLogLik = -1;
				probMaxNode.reset();

				for(ncomb = 0; ncomb < ncomblist; ncomb++) {

					nocache = 1;
					if(!pestim->m_pCacheMutex) {
						nocache = !getCachedProb(pcomblist[ncomb], d, nnode, idparset, d, &probMaxNode, &fMaxLogLik);
					}
					else {
						MUTEX_LOCK(pestim->m_pCacheMutex);
						nocache = !getCachedProb(pcomblist[ncomb], d, nnode, idparset, d, &probMaxNode, &fMaxLogLik);
						MUTEX_UNLOCK(pestim->m_pCacheMutex);
					}

					if(nocache) { 
						memcpy(idparset, pcomblist[ncomb], d*sizeof(int));
						// add pcomplist[j] parent set to nnode
						baseCatnet.setParents(nnode, idparset, d);
			     
						// add perturbation
						if(perturbations && psubsamples) {
							numsubsamples = 0;
							for(j = 0; j < numsamples; j++) {
								if(!perturbations[j * numnodes + nnode]) {
									memcpy(psubsamples + numsubsamples*numnodes, psamples + j*numnodes, numnodes*sizeof(int));
									numsubsamples++;
								}
							}
							fMaxLogLik = baseCatnet.setNodeSampleProb(nnode, psubsamples, numsubsamples);
						}
						else {
							fMaxLogLik = baseCatnet.setNodeSampleProb(nnode, psamples, numsamples);
						}

						pProbNode = baseCatnet.getNodeProb(nnode);
						if(pProbNode)
							probMaxNode = *pProbNode;
						
						if(pestim->m_pCacheMutex) {
							MUTEX_LOCK(pestim->m_pCacheMutex);
							setCachedProb(idparset, d, nnode, idparset, d,
									&probMaxNode, fMaxLogLik);
							MUTEX_UNLOCK(pestim->m_pCacheMutex);
						}
						else
							setCachedProb(idparset, d, nnode, idparset, d,
									&probMaxNode, fMaxLogLik);

					} /* if(!getCachedProb) */

					/* find nnode-complexity for pcomblist[ncomb] parent set*/
					nodecomplx = m_pNodeNumCats[nnode]-1;
					for(k = 0; k < d; k++)
						nodecomplx *= m_pNodeNumCats[idparset[k]];
					for(k = 0; k < m_nCatnets; k++) {
						if(!m_pCatnets[k])
							continue;
						pProbNode = m_pCatnets[k]->getNodeProb(nnode);
						complx = m_pCatnets[k]->complexity() - 
							m_pCatnets[k]->nodeComplexity(nnode) +
							nodecomplx;
						if(complx > maxComplexity || !pProbNode) 
							continue;
						tempLogLik = m_pCatnets[k]->loglik() - 
							pProbNode->loglik + fMaxLogLik;
						if(!pCurCatnetList[complx] && tempLogLik > -FLT_MAX) {
							pCurCatnetList[complx] = new CATNET<t_node, t_node_size, t_prob>;
							if (!pCurCatnetList[complx])
								return CATNET_ERR_MEM;
						}
						if(pCurCatnetList[complx] && 
							pCurCatnetList[complx]->loglik() < tempLogLik) {
								*pCurCatnetList[complx] = *m_pCatnets[k];
								pCurCatnetList[complx]->setParents(nnode, idparset, d);
								pCurCatnetList[complx]->setNodeProb(nnode, &probMaxNode);
						}
					} /* for k */

				} /* for ncomb */
	
				/* release combination set */
        			for(ncomb = 0; ncomb < ncomblist; ncomb++) {
          				if(pcomblist[ncomb])
            					CATNET_FREE(pcomblist[ncomb]);
          				pcomblist[ncomb] = NULL;
				} /* for ncomb */
        			CATNET_FREE(pcomblist);
        			pcomblist = 0;
				ncomblist = 0;

			} /* for d */

			} /* if(!bEqualCategories) */

			if(becho)
				Rprintf("\n");

			for(j = 0; j < m_nCatnets; j++) {
				if(m_pCatnets[j]) {
					if(pCurCatnetList[j]) {
						if(m_pCatnets[j]->loglik() < pCurCatnetList[j]->loglik()) {
							delete m_pCatnets[j];
							m_pCatnets[j] = pCurCatnetList[j];
							pCurCatnetList[j] = 0;
						}
						else {
							delete pCurCatnetList[j];
							pCurCatnetList[j] = 0;
						}
					}
				}
				else {
					m_pCatnets[j] = pCurCatnetList[j];
					pCurCatnetList[j] = 0;
				}
			}
		} // for(nnode = 0; nnode < numnodes; nnode++)

		CATNET_FREE(pCurCatnetList);
		CATNET_FREE(parset);
		CATNET_FREE(fixparset);
		CATNET_FREE(idparset);

		if(psubsamples)
			CATNET_FREE(psubsamples);

		for(j = 0; j < m_nCatnets; j++) {
			if(m_pCatnets[j]) {
				m_pCatnets[j] -> normalizeProbabilities(numsamples);
				m_pCatnets[j] -> setCategoryIndices(m_pNodeNumCats, m_pNodeCats);
			}
		} 

		if(m_pNodeCats) {
			for(i = 0; i < m_numNodes; i++) 
				if(m_pNodeCats[i])
					CATNET_FREE(m_pNodeCats[i]);
			CATNET_FREE(m_pNodeCats);
		}
		m_pNodeCats = 0;
		if(m_pNodeNumCats) 
			CATNET_FREE(m_pNodeNumCats);
		m_pNodeNumCats = 0;
		
		return CATNET_ERR_OK;
	}


};

#endif /* CATNET_SEARCH_H */
