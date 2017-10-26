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
 * catnet_class.h
 *
 *  Created on: Sep 18, 2009
 *      Author: nbalov
 *
 */

/* 
 * version 1.15.1  12dec2016
 */

#include "utils.h"
#include "problist.h"

#ifndef CATNET_CLASS_H_
#define CATNET_CLASS_H_

using namespace std;

template<class t_node, int t_node_size, class t_prob>
class CATNET {
protected:
	/* nodes are assumed ordered */
	int      m_numNodes;
	t_node **m_nodeNames;
	int      m_maxParents;
	int     *m_numParents;
	int    **m_parents;
	int      m_maxCategories;
	int     *m_numCategories;
	int    **m_catIndices;
	int      m_complexity;
	t_prob   m_loglik;
	PROB_LIST<t_prob> **m_pProbLists;

public:
	CATNET() {
		_reset();
	}

	CATNET(int nnodes, int maxpars, int maxcats = 2, const t_node **nodes = 0,
			const int * pnumpars = 0, const int **ppars = 0, const int *pcats = 0) {
		_reset();
		init(nnodes, maxpars, maxcats, nodes, pnumpars, ppars, pcats);
	}

	virtual ~CATNET() {
		_release();
	}

	CATNET<t_node, t_node_size, t_prob>& operator =(const CATNET<t_node,
			t_node_size, t_prob> &cnet) {
		init(cnet.m_numNodes, cnet.m_maxParents, cnet.m_maxCategories,
				(const t_node**)cnet.m_nodeNames, (const int*)cnet.m_numParents, (const int**)cnet.m_parents,
				(const int*)cnet.m_numCategories, (const PROB_LIST<t_prob> **)cnet.m_pProbLists);
		m_complexity = cnet.m_complexity;
		m_loglik = cnet.m_loglik;
		return *this;
	}

	virtual void _release() {
		int i;
		for (i = 0; i < m_numNodes; i++) {
			if (m_pProbLists && m_pProbLists[i]) {
				delete m_pProbLists[i];
				m_pProbLists[i] = 0;
			}
			if (m_parents && m_parents[i]) {
				CATNET_FREE(m_parents[i]);
				m_parents[i] = 0;
			}
			if (m_nodeNames && m_nodeNames[i]) {
				CATNET_FREE(m_nodeNames[i]);
				m_nodeNames[i] = 0;
			}
			if(m_catIndices && m_catIndices[i]) {
				CATNET_FREE(m_catIndices[i]);
				m_catIndices[i] = 0;
			}
		}
		if (m_numParents)
			CATNET_FREE(m_numParents);
		if (m_parents)
			CATNET_FREE(m_parents);
		if (m_numCategories)
			CATNET_FREE(m_numCategories);
		if (m_nodeNames)
			CATNET_FREE(m_nodeNames);
		if(m_catIndices)
			CATNET_FREE(m_catIndices);
		if (m_pProbLists)
			CATNET_FREE(m_pProbLists);

		_reset();
	}

	virtual void _reset() {
		m_numNodes = 0;
		m_maxParents = 0;
		m_maxCategories = 0;
		m_nodeNames = 0;
		m_numParents = 0;
		m_parents = 0;
		m_numCategories = 0;
		m_catIndices = 0;
		m_pProbLists = 0;
		m_complexity = 0;
		m_loglik = 0;
	}

	void init(int nnodes, int maxpars, int maxcats = 2, const t_node **nodes = 0,
			const int *pnumpars = 0, const int **ppars = 0, const int *pcats = 0, 
			const PROB_LIST<t_prob> **pprobs = 0) {

		_release();

		int i, j, *nodeparcats, nodenamelen;

		if (nnodes < 1 || maxpars < 0)
			return;

		m_numParents = (int*) CATNET_MALLOC(nnodes * sizeof(int));
		if (!m_numParents)
			return;
		m_parents = (int**) CATNET_MALLOC(nnodes * sizeof(int*));		
		if (!m_parents) {
			CATNET_FREE(m_numParents);
			m_numParents = 0;
			return;
		}
		m_numCategories = (int*) CATNET_MALLOC(nnodes * sizeof(int));
		if (!m_numCategories) {
			CATNET_FREE(m_numParents);
			m_numParents = 0;
			CATNET_FREE(m_parents);
			m_parents = 0;
			return;
		}

		m_numNodes = nnodes;
		m_maxParents = maxpars;
		m_maxCategories = maxcats;

		memset(m_numParents, 0, m_numNodes * sizeof(int));
		memset(m_parents,    0, m_numNodes * sizeof(int*));
		if (pcats) {
			memcpy(m_numCategories, pcats, m_numNodes * sizeof(int));
		}
		else {
			for (i = 0; i < m_numNodes; i++)
				m_numCategories[i] = m_maxCategories;
		}

		m_catIndices = (int**) CATNET_MALLOC(m_numNodes * sizeof(int*));
		if (m_catIndices)
			memset(m_catIndices, 0, m_numNodes * sizeof(int*));
		
		m_nodeNames = (t_node**) CATNET_MALLOC(m_numNodes * sizeof(t_node*));
		if (m_nodeNames) {
			if (nodes) {	
				for (i = 0; i < m_numNodes; i++) {
					m_nodeNames[i] = 0;
					if (!nodes[i])
						continue;
					nodenamelen = strlen(nodes[i]);
					m_nodeNames[i] = (t_node*) CATNET_MALLOC((nodenamelen+1) * sizeof(t_node));
					if (m_nodeNames[i] && nodes[i])
						strcpy(m_nodeNames[i], nodes[i]);
				}
			}
			else {
				memset(m_nodeNames, 0, m_numNodes * sizeof(t_node*));
			}
		}

		if (pnumpars && ppars) {
			memcpy(m_numParents, pnumpars, m_numNodes * sizeof(int));
			for (i = 0; i < m_numNodes; i++) {
				if (m_numParents[i] <= 0)
					continue;
				m_parents[i] = (int*) CATNET_MALLOC(m_numParents[i] * sizeof(int));
				if (m_parents[i]) {
					memset(m_parents[i], 0, m_numParents[i] * sizeof(int));
					if (ppars[i])
						memcpy(m_parents[i], ppars[i], m_numParents[i] * sizeof(int));
				}
				if (m_numParents[i] > m_maxParents)
					m_maxParents = m_numParents[i];
			}
		}
		
		m_pProbLists = (PROB_LIST<t_prob>**) CATNET_MALLOC(m_numNodes * sizeof(PROB_LIST<t_prob>*));
		if (!m_pProbLists) {
			CATNET_FREE(m_numParents);
			m_numParents = 0;
			CATNET_FREE(m_parents);
			m_parents = 0;
			CATNET_FREE(m_numCategories);
			m_numCategories = 0;
			if (m_catIndices) {
				CATNET_FREE(m_catIndices);
				m_catIndices = 0;
			}
			if (m_nodeNames) {
				CATNET_FREE(m_nodeNames);
				m_nodeNames = 0;
			}
			return;
		}
		memset(m_pProbLists, 0, m_numNodes * sizeof(PROB_LIST<t_prob>*));
		nodeparcats = (int*)CATNET_MALLOC((m_maxParents+1)*sizeof(int));
		if (nodeparcats) {			
			for (i = 0; i < m_numNodes; i++) {
				if (pprobs && pprobs[i]) {
					m_pProbLists[i] = new PROB_LIST<t_prob>;
					if (m_pProbLists[i])
						*m_pProbLists[i] = *pprobs[i];
				}
				else {
					for(j = 0; j < m_numParents[i]; j++)
						nodeparcats[j] = m_numCategories[m_parents[i][j]]; 
					m_pProbLists[i] = new PROB_LIST<t_prob>(m_numCategories[i], 								m_maxCategories, m_numParents[i], nodeparcats);
				}
			}
			CATNET_FREE(nodeparcats);
		}
	}

	void setNodeNames(char **pnames, const int *porder) {
		int i;
		const char *str;
		if (!porder || !pnames) 
			return;
		if(!m_nodeNames) {
			m_nodeNames = (t_node**) CATNET_MALLOC(m_numNodes * sizeof(t_node*));
			if (m_nodeNames)
				memset(m_nodeNames, 0, m_numNodes * sizeof(t_node*));
		}
		if(!m_nodeNames)
			return;
		for (i = 0; i < m_numNodes; i++) {
			m_nodeNames[i] = 0;
			if (porder[i] < 1 || porder[i] > m_numNodes)
				break;
			str = pnames[porder[i]-1];
			if(m_nodeNames[i])
				CATNET_FREE(m_nodeNames[i]);
			m_nodeNames[i] = (t_node*) CATNET_MALLOC((strlen(str)+1) * sizeof(char));
			if (m_nodeNames[i] && str)
				strcpy((char*)m_nodeNames[i], str);
		}
	}

	void setNodesOrder(const int *porder) {
		int i;
		char str[256];
		if (!porder) 
			return;
		if(!m_nodeNames) {
			m_nodeNames = (t_node**) CATNET_MALLOC(m_numNodes * sizeof(t_node*));
			if (m_nodeNames)
				memset(m_nodeNames, 0, m_numNodes * sizeof(t_node*));
		}
		if(!m_nodeNames)
			return;
		for (i = 0; i < m_numNodes; i++) {
			m_nodeNames[i] = 0;
			if (porder[i] < 1 || porder[i] > m_numNodes)
				break;
			sprintf(str, "N%d", (int)porder[i]);
			if(m_nodeNames[i])
				CATNET_FREE(m_nodeNames[i]);
			m_nodeNames[i] = (t_node*) CATNET_MALLOC((strlen(str)+1) * sizeof(char));
			if (m_nodeNames[i] && str)
				strcpy((char*)m_nodeNames[i], str);
		}
	}

	void setCategoryIndices(int *pNumCats, int **pCatIndices) {
		int i;
		// reset node names only
		if (!pNumCats || !pCatIndices || !m_numCategories) 
			return;
		if(!m_catIndices) {
			m_catIndices = (int**) CATNET_MALLOC(m_numNodes * sizeof(int*));
			if (m_catIndices)
				memset(m_catIndices, 0, m_numNodes * sizeof(int*));
		}
		if(!m_catIndices)
			return;
		for (i = 0; i < m_numNodes; i++) {
			if(m_numCategories[i] < 1 || pNumCats[i] < m_numCategories[i])
				continue;
			if(!m_catIndices[i])
				m_catIndices[i] = (int*) CATNET_MALLOC(m_numCategories[i] * sizeof(int));
			if (m_catIndices[i] && pCatIndices[i])
				memcpy(m_catIndices[i], pCatIndices[i], m_numCategories[i] * sizeof(int));
		}
	}

	void setCondProb(int node, t_prob *pcondprob, int condprobsize) {
		if (m_numNodes < 1 || node < 0 || node >= m_numNodes)
			return;
		if (!m_pProbLists) {
			m_pProbLists = (PROB_LIST<t_prob>**) CATNET_MALLOC(m_numNodes
					* sizeof(PROB_LIST<t_prob>*));
			if (m_pProbLists)
				memset(m_pProbLists, 0, m_numNodes * sizeof(PROB_LIST<t_prob>*));
		}
		if (!m_pProbLists)
			return;
		if (m_pProbLists[node])
			delete m_pProbLists[node];
		m_pProbLists[node] = 0;
		if (m_numParents[node] < 0 || m_numParents[node] > m_maxParents)
			return;
		int *parcats = (int*) CATNET_MALLOC(m_maxParents * sizeof(int));
		for (int i = 0; i < m_numParents[node]; i++) {
		        parcats[i] = m_numCategories[m_parents[node][i]];
		}
		m_pProbLists[node] = new PROB_LIST<t_prob> (m_numCategories[node],
				m_maxCategories, m_numParents[node], parcats, pcondprob,
				condprobsize);
		CATNET_FREE(parcats);
	}

	int numNodes() {
		return m_numNodes;
	}

	const t_node** nodeNames() {
		return (const t_node**)m_nodeNames;
	}

	int maxCategories() {
		return m_maxCategories;
	}

	const int* numCategories() {
		return (const int*)m_numCategories;
	}

	int numCategories(int node) {
		if (node < 0 || node >= m_numNodes)
			return 0;
		return m_numCategories[node];
	}

	int maxParents() {
		return m_maxParents;
	}

	const int** parents() {
		return (const int**)m_parents;
	}

	const int* numParents() {
		return (const int*)m_numParents;
	}

	int numParents(int node) {
		if (node < 0 || node >= m_numNodes)
			return 0;
		return m_numParents[node];
	}

	const int* getParents(int node) {
		if (node < 0 || node >= m_numNodes)
			return 0;
		return (const int*)m_parents[node];
	}

	int setParents(int node, int* parents, int numparents) {
		if (node < 0 || node >= m_numNodes)
			return 0;
		if (!m_numParents || !m_parents)
			return 0;
		if(m_numParents[node] != numparents) {
			m_numParents[node] = numparents;
			if(m_parents[node])
				CATNET_FREE(m_parents[node]);
			m_parents[node] = (int*) CATNET_MALLOC(m_numParents[node] * sizeof(int));
		}
		if (m_parents[node] && parents)
			memcpy(m_parents[node], parents, m_numParents[node] * sizeof(int));

		if (!m_pProbLists) {
			m_pProbLists = (PROB_LIST<t_prob>**) CATNET_MALLOC(m_numNodes * sizeof(PROB_LIST<t_prob>*));
			if (!m_pProbLists)
				return 0;
			memset(m_pProbLists, 0, m_numNodes * sizeof(PROB_LIST<t_prob>*));
		}
		else {
			if(m_pProbLists[node])
				delete m_pProbLists[node];
		}
		int *parcats = (int*)CATNET_MALLOC(m_numParents[node]*sizeof(int));
		if (!parcats)
			return 0;
		for(int i = 0; i < m_numParents[node]; i++) 
			parcats[i] = m_numCategories[parents[i]];	
		m_pProbLists[node] = new PROB_LIST<t_prob>(m_numCategories[node], m_maxCategories, m_numParents[node], parcats);
		CATNET_FREE(parcats);
		             
		if(m_maxParents < numparents)
			m_maxParents = numparents;

		/* need to be calculated next time */
		m_complexity = 0;
		m_loglik = 0;

		return m_numParents[node];
	}

	int complexity() {
		if(m_complexity < m_numNodes)
			return findComplexity();
		return m_complexity;
	}

	int findComplexity() {
		int i, j, c;
		m_complexity = 0;
		for (i = 0; i < m_numNodes; i++) {
			if (!m_parents || !m_parents[i]) {
				m_complexity += (m_numCategories[i]-1);
				continue;
			}
			c = 1;
			for (j = 0; j < m_numParents[i]; j++)
				c *= m_numCategories[m_parents[i][j]];
			m_complexity += c*(m_numCategories[i]-1);
		}
		return m_complexity;
	}

	int nodeComplexity(int nnode) {
		int j, c;
		if(nnode < 0 || nnode >= m_numNodes)
			return(0);
		if (!m_parents || !m_parents[nnode])
			return(m_numCategories[nnode]-1);
		c = (m_numCategories[nnode]-1);
		for (j = 0; j < m_numParents[nnode]; j++)
			c *= m_numCategories[m_parents[nnode][j]];
		return c;
	}

	t_prob loglik() {
		if(m_loglik == 0)
			return findLoglik();
		return m_loglik;
	}

	t_prob findLoglik() {
		int i;
		if(!m_pProbLists)
			return -FLT_MAX;
		m_loglik = 0;
		for (i = 0; i < m_numNodes; i++) {
			if(m_pProbLists[i])
				m_loglik += (m_pProbLists[i]->loglik + m_pProbLists[i]->priorlik);
		}
		return m_loglik;
	}
	
	const PROB_LIST<t_prob>** probLists() {
		return (const PROB_LIST<t_prob>**)m_pProbLists;
	}

	PROB_LIST<t_prob>* getNodeProb(int nnode) {
		if (!m_pProbLists || nnode < 0 || nnode >= m_numNodes)
			return 0;
		return m_pProbLists[nnode];
	}

	const PROB_LIST<t_prob>* setNodeProb(int nnode, const PROB_LIST<t_prob> *pprob) {
		if(!m_pProbLists || nnode < 0 || nnode >= m_numNodes)
			return 0;
		if(!m_pProbLists[nnode])
			m_pProbLists[nnode] = new PROB_LIST<t_prob>;
		if (m_pProbLists[nnode] && pprob) {
			*m_pProbLists[nnode] = *pprob;
		}
		return m_pProbLists[nnode];
	}

	void setNodePriorProb(int nnode, t_prob pprior) {
		if(!m_pProbLists || nnode < 0 || nnode >= m_numNodes || !m_pProbLists[nnode])
			return;
		m_pProbLists[nnode]->priorlik = pprior;
		m_loglik = 0;
	}

	t_prob nodeSampleLoglik(int nnode, int *pnodepars, int nodepars,
			int *psamples, int nsamples) {
		int j, ipar;
		/* psamples have categories in the range [1, m_maxCategories] */
		t_prob *pnodeprob, floglik = 0;
		int *pnodesample, samp, ncount;

		if(!psamples || nsamples < 1)
			return 0;

		pnodesample = 0;
		if (m_maxParents > 0) {
			pnodesample = (int*) CATNET_MALLOC(m_maxParents * sizeof(int));
		}
		pnodepars = m_parents[nnode];
		ncount = 0;
		for (j = 0; j < nsamples; j++) {
			if (pnodesample) {
				for (ipar = 0; ipar < nodepars; ipar++) {
					pnodesample[ipar] = psamples[j * m_numNodes + pnodepars[ipar]];
				}
			}
			// for pnodesample = 0, pnodeprob = m_pProbLists[nnode]->pProbs
			pnodeprob = m_pProbLists[nnode]->find_slot(0, pnodesample, 0);
			samp = psamples[j * m_numNodes + nnode];
			if (samp >= 0 && samp < m_numCategories[nnode]) {
				if(pnodeprob[samp] > 0)
					floglik += (t_prob)log((double)pnodeprob[samp]);
				else {
					floglik = (t_prob)-FLT_MAX;
					break;
				}
				ncount++;
			}
		}
		if(ncount > 1 && floglik > (t_prob)-FLT_MAX)
			floglik /= (t_prob)ncount;
		CATNET_FREE(pnodesample);
		return(floglik);
	}

	t_prob sampleLoglik(int *psamples, int nsamples) {
		int i, j, ipar;
		/* psamples have categories in the range [1, m_maxCategories] */
		t_prob *pnodeprob, nodeloglik;
		int nodepars, ncount;
		int *pnodepars, *pnodesample, samp;

		if(!psamples || nsamples < 1)
			return 0;

		pnodesample = 0;
		if (m_maxParents > 0) {
			pnodesample = (int*) CATNET_MALLOC(m_maxParents * sizeof(int));
		}
		m_loglik = 0;
		for (i = 0; i < m_numNodes; i++) {
			if(!m_pProbLists[i])
				continue;
			pnodepars  = m_parents[i];
			nodepars   = m_numParents[i];
			ncount     = 0;
			nodeloglik = 0;
			for (j = 0; j < nsamples; j++) {
				if (pnodesample) {
					for (ipar = 0; ipar < nodepars; ipar++) {
						pnodesample[ipar] = psamples[j * m_numNodes + pnodepars[ipar]];
					}
				}
				// for pnodesample = 0, pnodeprob = m_pProbLists[nnode]->pProbs
				pnodeprob = m_pProbLists[i]->find_slot(0, pnodesample, 0);
				samp = psamples[j * m_numNodes + i];
				if (pnodeprob && samp >= 0 && samp < m_numCategories[i]) {
					if(pnodeprob[samp] > 0)
						nodeloglik += (t_prob)log((double)pnodeprob[samp]);
					else {
						nodeloglik = (t_prob)-FLT_MAX;
						break;
					}
					ncount++;
				}
			}
			if(ncount > 1 && nodeloglik > (t_prob)-FLT_MAX)
				nodeloglik /= (t_prob)ncount;
			if(nodeloglik == -FLT_MAX) {
				m_loglik = -FLT_MAX;
				break;
			}
			else
				m_loglik += nodeloglik;
		}
		CATNET_FREE(pnodesample);
		return m_loglik;
	}

	t_prob* sampleLoglikVector(int *psamples, int nsamples, int *pert=0) {
		int i, j, ipar;
		/* psamples have categories in the range [1, m_maxCategories] */
		t_prob *pnodeprob, *ploglik;
		int nodepars, ncount;
		int *pnodepars, *pnodesample, samp;

		if(!psamples || nsamples < 1)
			return 0;

		ploglik = (t_prob*) CATNET_MALLOC(m_numNodes * sizeof(t_prob));
		if (!ploglik)
			return 0;
		memset(ploglik, 0, m_numNodes * sizeof(t_prob));

		pnodesample = 0;
		if (m_maxParents > 0) {
			pnodesample = (int*) CATNET_MALLOC(m_maxParents * sizeof(int));
		}
		for (i = 0; i < m_numNodes; i++) {
			if(!m_pProbLists[i])
				continue;
			pnodepars = m_parents[i];
			nodepars = m_numParents[i];

			ncount = 0;
			if(pert){
				for (j = 0; j < nsamples; j++) {
					// check for perturbation
					if(pert[j * m_numNodes + i])
						continue; 
					if (pnodesample) {
						for (ipar = 0; ipar < nodepars; ipar++) {
							pnodesample[ipar] = psamples[j * m_numNodes + pnodepars[ipar]];
						}
					}
					// for pnodesample = 0, pnodeprob = m_pProbLists[nnode]->pProbs
					pnodeprob = m_pProbLists[i]->find_slot(0, pnodesample, 0);
					samp = psamples[j * m_numNodes + i];
					if (pnodeprob && samp >= 0 && samp < m_numCategories[i]) {
						if(pnodeprob[samp] > 0)
							ploglik[i] += (t_prob)log((double)pnodeprob[samp]);
						else {
							ploglik[i] = (t_prob)-FLT_MAX;
							break;
						}
						ncount++;
					}
				}
			}
			else {
				for (j = 0; j < nsamples; j++) {
					if (pnodesample) {
						for (ipar = 0; ipar < nodepars; ipar++) {
							pnodesample[ipar] = psamples[j * m_numNodes + pnodepars[ipar]];
						}
					}
					pnodeprob = m_pProbLists[i]->find_slot(0, pnodesample, 0);
					samp = psamples[j * m_numNodes + i];
					if (pnodeprob && samp >= 0 && samp < m_numCategories[i]) {
						if(pnodeprob[samp] > 0)
							ploglik[i] += (t_prob)log((double)pnodeprob[samp]);
						else {
							ploglik[i] = (t_prob)-FLT_MAX;
							break;
						}
						ncount++;
					}
				}
			}
			if(ncount > 1 && ploglik[i] > -FLT_MAX)
				ploglik[i] /= (t_prob)ncount;
		}
		CATNET_FREE(pnodesample);
		return ploglik;
	}

	t_prob* bySampleLoglikVector(int *psamples, int nsamples, int *pert=0) {
		int i, j, ipar;
		/* psamples have categories in the range [1, m_maxCategories] */
		t_prob *pnodeprob, *ploglik;
		int nodepars, *pncount;
		int *pnodepars, *pnodesample, samp;

		if(!psamples || nsamples < 1)
			return 0;

		ploglik = (t_prob*) CATNET_MALLOC(nsamples * sizeof(t_prob));
		if (!ploglik)
			return 0;
		memset(ploglik, 0, nsamples * sizeof(t_prob));

		pncount = (int*) CATNET_MALLOC(nsamples * sizeof(int));
		if (!pncount) {
			CATNET_FREE(ploglik);
			return 0;
		}

		pnodesample = 0;
		if (m_maxParents > 0) {
			pnodesample = (int*) CATNET_MALLOC(m_maxParents * sizeof(int));
		}
		for (i = 0; i < m_numNodes; i++) {
			if(!m_pProbLists[i])
				continue;
			pnodepars = m_parents[i];
			nodepars = m_numParents[i];

			if(pert){
				for (j = 0; j < nsamples; j++) {
					// check for perturbation
					if(pert[j * m_numNodes + i])
						continue; 
					if (pnodesample) {
						for (ipar = 0; ipar < nodepars; ipar++) {
							pnodesample[ipar] = psamples[j * m_numNodes + pnodepars[ipar]];
						}
					}
					// for pnodesample = 0, pnodeprob = m_pProbLists[nnode]->pProbs
					pnodeprob = m_pProbLists[i]->find_slot(0, pnodesample, 0);
					samp = psamples[j * m_numNodes + i];
					if (pnodeprob && samp >= 0 && samp < m_numCategories[i]) {
						if(pnodeprob[samp] > 0)
							ploglik[j] += (t_prob)log((double)pnodeprob[samp]);
						else {
							ploglik[j] = (t_prob)-FLT_MAX;
							break;
						}
						pncount[j]++;
					}
				}
			}
			else {
				for (j = 0; j < nsamples; j++) {
					for (ipar = 0; ipar < nodepars; ipar++) {
						pnodesample[ipar] = psamples[j * m_numNodes + pnodepars[ipar]];
					}
					pnodeprob = m_pProbLists[i]->find_slot(0, pnodesample, 0);
					samp = psamples[j * m_numNodes + i];
					if (pnodeprob && samp >= 0 && samp < m_numCategories[i]) {
						if(pnodeprob[samp] > 0)
							ploglik[j] += (t_prob)log((double)pnodeprob[samp]);
						else {
							ploglik[j] = (t_prob)-FLT_MAX;
							break;
						}
						pncount[j]++;
					}
				}
			}
		}

		CATNET_FREE(pnodesample);
		CATNET_FREE(pncount);
		return ploglik;
	}

	t_prob sampleNodeLoglik(int nnode, int *psamples, int nsamples) {
		int j, ipar;
		/* psamples have categories in the range [1, m_maxCategories] */
		t_prob *pnodeprob, loglik;
		int nodepars, ncount;
		int *pnodepars, *pnodesample=0, samp;

		if(!psamples || nsamples < 1 || nnode < 0 || nnode >= m_numNodes)
			return 0;

		if(!m_pProbLists || !m_pProbLists[nnode])
			return 0;

		loglik = 0;

		pnodesample = 0;
		if (m_maxParents > 0) {
			pnodesample = (int*) CATNET_MALLOC(m_maxParents * sizeof(int));	
		}
		pnodepars = m_parents[nnode];
		nodepars  = m_numParents[nnode];
	
		ncount = 0;	
		for (j = 0; j < nsamples; j++) {
			if (pnodesample) {
				for (ipar = 0; ipar < nodepars; ipar++) {
					pnodesample[ipar] = psamples[j * m_numNodes + pnodepars[ipar]];
				}
			}
			// for pnodesample = 0, pnodeprob = m_pProbLists[nnode]->pProbs
			pnodeprob = m_pProbLists[nnode]->find_slot(0, pnodesample, 0);
			samp = psamples[j * m_numNodes + nnode];
			if (pnodeprob && samp >= 0 && samp < m_numCategories[nnode]) {
				if(pnodeprob[samp] > 0)
					loglik += (t_prob)log((double)pnodeprob[samp]);
				else {
					loglik = (t_prob)-FLT_MAX;
					break;
				}
				ncount++;
			}
		}
		CATNET_FREE(pnodesample);

		if(ncount > 1 && loglik > (t_prob)-FLT_MAX)
			loglik /= (t_prob)ncount;

		return(loglik);
	}

	// sets sample conditional probability and returns its log-likelihood
	t_prob setNodeSampleProb(int nnode, int *psamples, int nsamples, int bNormalize = 0) {
		int i, j;
		/* psamples have categories in the range [1, m_maxCategories] */
		t_prob *pnodeprob, floglik;
		int *pnodesample, *pnodepars, samp, ncount;

		if (!m_pProbLists || !psamples || nsamples < 1) {
			return (t_prob)-FLT_MAX;
		}

		if(!m_pProbLists[nnode]) {
			// error
			return (t_prob)-FLT_MAX;
		}
		else {
			m_pProbLists[nnode]->set_zero();
		}

		pnodesample = 0;
		if (m_maxParents > 0) {
			pnodesample = (int*) CATNET_MALLOC(m_maxParents * sizeof(int));
		}
		pnodepars = m_parents[nnode];
		ncount = 0;
		for (j = 0; j < nsamples; j++) {
			if (pnodesample) {
				for (i = 0; i < m_numParents[nnode]; i++) {
					if (pnodepars[i] < 0 || pnodepars[i] >= m_numNodes)
						break;
					pnodesample[i] = psamples[j * m_numNodes + pnodepars[i]];
				}
			}
			// for pnodesample = 0, pnodeprob = m_pProbLists[nnode]->pProbs
			pnodeprob = m_pProbLists[nnode]->find_slot(0, pnodesample, 0);
			samp = psamples[j * m_numNodes + nnode];
			if (pnodeprob && samp >= 0 && samp < m_numCategories[nnode]) {
				pnodeprob[samp]++;
				ncount++;
			}
		}

		CATNET_FREE(pnodesample);

		/* keep sample sizes */
		m_pProbLists[nnode]->sampleSize = ncount;

		/* at this point m_pProbLists[nnode] has counts not probabilities 
		  find m_pProbLists[nnode]->loglik */
		m_pProbLists[nnode]->loglikelihood();
		if(ncount > 1) {
			m_pProbLists[nnode]->loglik   /= (t_prob)ncount;
			m_pProbLists[nnode]->priorlik /= (t_prob)ncount;
		}
		floglik = m_pProbLists[nnode]->loglik + m_pProbLists[nnode]->priorlik;
		if(bNormalize)
			m_pProbLists[nnode] -> normalize();

		/* need to be calculated next call time */
		m_loglik = 0;

		return(floglik);
	}

	int normalizeProbabilities() {
		int i;
		if(!m_pProbLists)
			return -1;
		for (i = 0; i < m_numNodes; i++) {
			if(m_pProbLists[i]) {
				m_pProbLists[i] -> normalize();
			}
		}
		return 0;
	}

	int *findParentPool(int nnode, int &poolsize) {

		int ipar, par, i, j, bfound, parpoolsize, *parpool, *ppool, *paux;

		poolsize = 0;
		if (nnode < 0 || nnode >= m_numNodes || !m_parents || !m_parents[nnode] || !m_numParents[nnode])
			return 0;

		ppool = 0;
		for (ipar = 0; ipar < m_numParents[nnode]; ipar++) {
			par = m_parents[nnode][ipar];
			parpool = findParentPool(par, parpoolsize);
			if (parpool && parpoolsize > 0) {
				i = 0;
				while (i < parpoolsize) {
					bfound = 0;
					for (j = 0; j < poolsize; j++) {
						if (parpool[i] == ppool[j]) {
							bfound = 1;
							break;
						}
					}
					if (!bfound) {
						i++;
						continue;
					}
					for (j = i + 1; j < parpoolsize; j++)
						parpool[j - 1] = parpool[j];
					parpoolsize--;
				}
			}

			paux = (int*) CATNET_MALLOC((poolsize + parpoolsize + 1) * sizeof(int));
			if (!paux) {
				if(parpool)
					CATNET_FREE(parpool);
				return 0;
			}
			if (ppool && poolsize > 0)
				memcpy(paux, ppool, poolsize * sizeof(int));
			if (parpool && parpoolsize > 0)
				memcpy(paux + poolsize, parpool, parpoolsize * sizeof(int));
			// release that
			if(parpool)
				CATNET_FREE(parpool);
			parpool = 0;

			// check par
			bfound = 0;
			for (j = 0; j < poolsize + parpoolsize; j++) {
				if (paux[j] == par) {
					bfound = 1;
					break;
				}
			}
			if (bfound) 
				poolsize += parpoolsize;
			else {
				paux[poolsize + parpoolsize] = par;
				poolsize += (parpoolsize + 1);
			}

			if (ppool)
				CATNET_FREE(ppool);
			ppool = paux;

		}
		return ppool;
	}

	t_prob *findJointProb(int nnode, int &jointprobsize) {
		t_prob *jointprob, *parprob;
		int i, ii, i0, ic, cc, ipar, j, par, ipool, poolnode;
		int *parpool, *paux, parpoolsize, *blocksize, *paridx, *pcats;
		PROB_LIST<t_prob> *probnode = 0;

		if (nnode < 0 || nnode >= m_numNodes || !m_pProbLists)
			return 0;

		parpool = findParentPool(nnode, parpoolsize);
		// add nnode to the parents list
		paux = (int*) CATNET_MALLOC((parpoolsize + 1) * sizeof(int));
		if (!paux) 
			return 0;
		if (parpool) {
			if (parpoolsize > 0) 
				memcpy(paux, parpool, parpoolsize * sizeof(int));	
			CATNET_FREE(parpool);
		}
		paux[parpoolsize] = nnode;
		parpoolsize++;
		parpool = paux;

		pcats = 0;
		if (m_maxParents > 0) {
			pcats = (int*) CATNET_MALLOC(m_maxParents * sizeof(int));
			if (!pcats)
				return 0;
		}
		paridx = (int*) CATNET_MALLOC(parpoolsize * sizeof(int));
		if (!paridx) {
			CATNET_FREE(parpool);
			CATNET_FREE(pcats);
			return 0;
		}
		blocksize = (int*) CATNET_MALLOC(parpoolsize * sizeof(int));
		if (!blocksize) {
			CATNET_FREE(parpool);
			CATNET_FREE(pcats);
			CATNET_FREE(paridx);
			return 0;
		}

		blocksize[parpoolsize - 1] = 1;
		for (i = parpoolsize - 2; i >= 0; i--) {
			blocksize[i] = blocksize[i + 1] * m_numCategories[parpool[i + 1]];
		}

		jointprobsize = 1;
		for (i = 0; i < parpoolsize; i++) {
			jointprobsize *= m_numCategories[parpool[i]];
		}

		jointprob = (t_prob*) CATNET_MALLOC(jointprobsize * sizeof(t_prob));
		if (!jointprob) {
			CATNET_FREE(pcats);
			CATNET_FREE(paridx);
			CATNET_FREE(parpool);
			CATNET_FREE(blocksize);
			return 0;
		}
		for (i = 0; i < jointprobsize; i++)
			jointprob[i] = 1.0;

		for (ipool = 0; ipool < parpoolsize; ipool++) {
			poolnode = parpool[ipool];
			probnode = m_pProbLists[poolnode];

			if (m_numParents[poolnode] == 0) {
				for (ii = 0; ii < jointprobsize; ii += (blocksize[ipool]
						* m_numCategories[poolnode])) {
					for (ic = 0; ic < m_numCategories[poolnode]; ic++) {
						i0 = ic * blocksize[ipool];
						for (i = 0; i < blocksize[ipool]; i++) {
							jointprob[ii + i0 + i] *= probnode->pProbs[ic];
						}
					}
				}
				continue;
			}

			//cout << "Parents: ";
			if (paridx && parpoolsize > 0)
				memset(paridx, 0, parpoolsize * sizeof(int));
			for (ipar = 0; ipar < m_numParents[poolnode]; ipar++) {
				par = m_parents[poolnode][ipar];
				for (i = 0; i < ipool; i++)
					if (par == parpool[i])
						break;
				if (i >= ipool) {
					// bad
					break;
				}
				paridx[i] = ipar + 1;
			}
			if (pcats && m_maxParents > 0) {
				memset(pcats, 0, m_maxParents * sizeof(int));
			}
			for (ii = 0; ii < jointprobsize; ii += (blocksize[ipool] * m_numCategories[poolnode])) {

				for (j = 0; j < ipool; j++) {
					if (paridx[j] > 0) {
						ic = (int) (ii / blocksize[j]);
						cc = m_numCategories[m_parents[poolnode][paridx[j] - 1]];
						ic -= cc * (int) (ic / cc);
						pcats[paridx[j] - 1] = ic;
						if(pcats[paridx[j] - 1] >= m_numCategories[parpool[j]])
							pcats[paridx[j] - 1] = m_numCategories[parpool[j]] - 1;
					}
				}
				parprob = probnode->find_slot(0, pcats, 0);
				if (!parprob)
					continue;
				for (ic = 0; ic < m_numCategories[poolnode]; ic++) {
					i0 = ic * blocksize[ipool];
					for (i = 0; i < blocksize[ipool]; i++) {
						jointprob[ii + i0 + i] *= parprob[ic];
					}
				}
			}
		}

		CATNET_FREE(pcats);
		CATNET_FREE(paridx);
		CATNET_FREE(parpool);
		CATNET_FREE(blocksize);

		return jointprob;
	}

	t_prob *marginal_prob(int nnode) {

		t_prob *jointprob, *margprob;
		int jointprobsize = 0, nodecats, ic, k;
		nodecats = m_numCategories[nnode];
		margprob = (t_prob*) CATNET_MALLOC(nodecats * sizeof(t_prob));
		if (!margprob)
			return 0;
		jointprob = findJointProb(nnode, jointprobsize);
		if (!jointprob) {
			CATNET_FREE(margprob);
			return 0;
		}
		for (ic = 0; ic < nodecats; ic++) {
			k = ic;
			margprob[ic] = 0;
			while (k < jointprobsize) {
				margprob[ic] += jointprob[k];
				k += nodecats;
			}
		}
		CATNET_FREE(jointprob);

		return margprob;
	}

	int *getOrder() {
		int *porder = 0, *pfound = 0, bhaspar, i, j, k;
		if(m_numNodes < 1 || !m_numParents || !m_parents)
			return 0;
		porder = (int*)CATNET_MALLOC(m_numNodes * sizeof(int));
		if (!porder)
			return 0;
		pfound = (int*)CATNET_MALLOC(m_numNodes * sizeof(int));
		if (!pfound) {
			CATNET_FREE(porder);
			return 0;
		}
		memset(pfound, 0, m_numNodes * sizeof(int));
		for(i = 0; i < m_numNodes; i++) {
			for(j = 0; j < m_numNodes; j++) {
				if(pfound[j])
					continue;
				if(m_numParents[j] <= 0)
					break;
				bhaspar = 0;
				for(k = 0; k < m_numParents[j]; k++) {
					if(pfound[m_parents[j][k]])
						continue;
					bhaspar = 1;
					break;
				}
				if(!bhaspar)
					break;
			}
			if(j >= m_numNodes || pfound[j]) {
				CATNET_FREE(pfound);
				CATNET_FREE(porder);
				return 0;
			}
			porder[i] = j;
			pfound[j] = 1;
		}
		CATNET_FREE(pfound);
		return porder;
	}

};

#endif /* CATNET_CLASS_H_ */
