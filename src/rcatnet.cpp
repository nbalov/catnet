/*
 *  catnet : categorical Bayesian network inference
 *  Copyright (C) 2009--2017  Nikolay Balov
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
 * version 1.15.7  09mar2020
 */

#include "utils.h"
#include "catnet_class.h"
#include "rcatnet.h"

using namespace std;

RCatnet::RCatnet() {
}

RCatnet::RCatnet(SEXP cnet) {

	double *pvec;
	int nnode, i, nvec, *pn;
	char const *pstr;

	SEXP rname, rnodes, rparents, rcatlist, rproblist, pf, nodepars, nodeproblist, pint, rnodeprob;

	if(!isS4(cnet))
		return;

	rname     = GET_SLOT(cnet, install("objectName"));
	rnodes    = GET_SLOT(cnet, install("nodes"));
	rparents  = GET_SLOT(cnet, install("parents"));
	rcatlist  = GET_SLOT(cnet, install("categories"));
	rproblist = GET_SLOT(cnet, install("probabilities"));

	if(rnodes == R_NilValue || rparents == R_NilValue || rcatlist == R_NilValue || rproblist == R_NilValue) {
		return;
	}

	PROTECT(rname = AS_CHARACTER(rname));
	PROTECT(rnodes = AS_LIST(rnodes));
	PROTECT(rparents = AS_LIST(rparents));
	PROTECT(rcatlist = AS_LIST(rcatlist));
	PROTECT(rproblist = AS_LIST(rproblist));

	pint = GET_SLOT(cnet, install("numnodes"));
	m_numNodes = INTEGER_POINTER(pint)[0];

	pint = GET_SLOT(cnet, install("maxParents"));
	m_maxParents = INTEGER_POINTER(pint)[0];
	pint = GET_SLOT(cnet, install("maxCategories"));
	m_maxCategories = INTEGER_POINTER(pint)[0];

	pint = GET_SLOT(cnet, install("complexity"));
	m_complexity = INTEGER_POINTER(pint)[0];
	pint = GET_SLOT(cnet, install("likelihood"));
	m_loglik = NUMERIC_POINTER(pint)[0];

	if (length(rproblist) != m_numNodes) {
		UNPROTECT(5);
		warning("length(rproblist) != m_numNodes");
		return;
	}
	m_nodeNames = (char**) CATNET_MALLOC(m_numNodes * sizeof(char*));
	m_numParents = (int*) CATNET_MALLOC(m_numNodes * sizeof(int));
	m_parents = (int**) CATNET_MALLOC(m_numNodes * sizeof(int*));
	m_numCategories = (int*) CATNET_MALLOC(m_numNodes * sizeof(int));
	m_pProbLists = (PROB_LIST<double>**) CATNET_MALLOC(m_numNodes * sizeof(PROB_LIST<double>*));

	if (!m_nodeNames || !m_numParents || !m_parents || !m_numCategories || !m_pProbLists) {
		if (m_nodeNames)
			CATNET_FREE(m_nodeNames);
		m_nodeNames = 0;
		if (m_numParents)
			CATNET_FREE(m_numParents);
		m_numParents = 0;
		if (m_parents)
			CATNET_FREE(m_parents);
		m_parents = 0;
		if (m_numCategories)
			CATNET_FREE(m_numCategories);
		m_numCategories = 0;
		if (m_pProbLists)
			CATNET_FREE(m_pProbLists);
		m_pProbLists = 0;
		UNPROTECT(5);
		return;
	}

	memset(m_nodeNames,     0, m_numNodes * sizeof(char*));
	memset(m_numParents,    0, m_numNodes * sizeof(int));
	memset(m_parents,       0, m_numNodes * sizeof(int*));
	memset(m_numCategories, 0, m_numNodes * sizeof(int));

	for (i = 0; i < m_numNodes; i++)
		m_pProbLists[i] = 0;
	
	for(nnode = 0; nnode < m_numNodes; nnode++) {
		pf = VECTOR_ELT(rnodes, nnode);
		m_nodeNames[nnode] = 0;
		if(IS_VECTOR(pf)) {
			pstr = CHAR(asChar(pf));
			if(strlen(pstr) >= MAX_NODE_NAME) {
				m_nodeNames[nnode] = (char*) CATNET_MALLOC((strlen(pstr)+1) * sizeof(char));
				if (m_nodeNames[nnode] && pstr) 
					strcpy(m_nodeNames[nnode], pstr);
			}
			else {
				m_nodeNames[nnode] = (char*) CATNET_MALLOC(MAX_NODE_NAME * sizeof(char));
				if (m_nodeNames[nnode])
					sprintf(m_nodeNames[nnode], "N%d", nnode);
			}
		}

		pf = VECTOR_ELT(rparents, nnode);
		m_numParents[nnode] = 0;
		m_parents[nnode] = 0;
		if (IS_VECTOR(pf)) {
			m_numParents[nnode] = length(pf);
			pn = INTEGER_POINTER(pf);
			m_parents[nnode] = (int*) CATNET_MALLOC(m_numParents[nnode] * sizeof(int));
			for(i = 0; i < m_numParents[nnode]; i++) {
				m_parents[nnode][i] = pn[i] - 1;
			}
		}

		pf = VECTOR_ELT(rcatlist, nnode);
		m_numCategories[nnode] = length(pf);
	}

	// get probabilities
	if(!strcmp(CHAR(asChar(rname)), "catNetworkC")) {
		for (nnode = 0; nnode < m_numNodes; nnode++) {
			rnodeprob = VECTOR_ELT(rproblist, nnode);
			setCondProb(nnode, NUMERIC_POINTER(rnodeprob), length(rnodeprob));   
		}
	}
	else {
		for (nnode = 0; nnode < m_numNodes; nnode++) {
			nodepars = VECTOR_ELT(rparents, nnode);
			nodeproblist = VECTOR_ELT(rproblist, nnode);
			pvec = 0;
			nvec = 0;
			gen_prob_vector(nnode, nodepars, 0, rcatlist, nodeproblist, pvec, nvec);
			setCondProb(nnode, pvec, nvec);
			CATNET_FREE(pvec);
		}
	}

	UNPROTECT(5);
}

SEXP RCatnet::genRcatnet(const char * objectName = (const char*)"catNetwork") {

	char str[256];
	int node, i, *pslotcats, *pn;
	double *pf, floglik;
	SEXP plist, ppars, pcats, pnodeprobs, strnames, pint, pstr;

	if(strcmp(objectName, "catNetwork") && strcmp(objectName, "catNetworkC") ) 
		return R_NilValue;

	SEXP cnet = PROTECT(NEW_OBJECT(MAKE_CLASS(objectName)));

	PROTECT(plist = allocVector(STRSXP, 1));
	SET_STRING_ELT(plist, 0, mkChar(objectName));
	SET_SLOT(cnet, install("objectName"), plist);
	UNPROTECT(1);

	PROTECT(pint = NEW_INTEGER(1));
	INTEGER_POINTER(pint)[0] = m_numNodes;
	SET_SLOT(cnet, install("numnodes"), pint);
	UNPROTECT(1);

	PROTECT(strnames = allocVector(STRSXP, m_numNodes));
	for(node = 0; node < m_numNodes; node++) {
		if(m_nodeNames && m_nodeNames[node]) {
			SET_STRING_ELT(strnames, node, mkChar(m_nodeNames[node]));
		}
		else {
			sprintf(str, "N%d", node+1);
			SET_STRING_ELT(strnames, node, mkChar(str));
		}
	}
	SET_SLOT(cnet, install("nodes"), strnames);
 	UNPROTECT(1);

	PROTECT(pint = NEW_INTEGER(1));
	INTEGER_POINTER(pint)[0] = m_maxParents;
	SET_SLOT(cnet, install("maxParents"), pint);
	UNPROTECT(1);

	PROTECT(plist = allocVector(VECSXP, m_numNodes));
	for(node = 0; node < m_numNodes; node++) {
		if(m_numParents[node] <= 0) {
			SET_VECTOR_ELT(plist, node, R_NilValue);
			continue;
		}
		PROTECT(ppars = NEW_INTEGER(m_numParents[node]));
		pn = INTEGER_POINTER(ppars);
		for(i = 0; i < m_numParents[node]; i++)
			// remember to increase the index by 1
			pn[i] = m_parents[node][i] + 1;
		SET_VECTOR_ELT(plist, node, ppars);
		UNPROTECT(1);
	}
	SET_SLOT(cnet, install("parents"), plist);
	UNPROTECT(1);

	PROTECT(pint = NEW_INTEGER(1));
	INTEGER_POINTER(pint)[0] = m_maxCategories;
	SET_SLOT(cnet, install("maxCategories"), pint);
	UNPROTECT(1);

	PROTECT(plist = allocVector(VECSXP, m_numNodes));
	for(node = 0; node < m_numNodes; node++) {
		if(m_numCategories[node] > m_maxCategories)
			break;
		PROTECT(pcats = allocVector(STRSXP, m_numCategories[node]));
		for(i = 0; i < m_numCategories[node]; i++) {
			if(m_catIndices && m_catIndices[node])
				sprintf(str, "%d", m_catIndices[node][i]);
			else
				sprintf(str, "C%d", i+1);
			SET_STRING_ELT(pcats, i, mkChar(str));
		}
		SET_VECTOR_ELT(plist, node, pcats);
		UNPROTECT(1);
	}
	SET_SLOT(cnet, install("categories"), plist);
	UNPROTECT(1); /* plist */

	pslotcats = (int*)CATNET_MALLOC(m_maxParents*sizeof(int));
	PROTECT(plist = allocVector(VECSXP, m_numNodes));
	if(!strcmp(objectName, "catNetworkC")) {	
		for(node = 0; node < m_numNodes; node++) {
			if(!m_pProbLists[node])
				continue;
			PROTECT(pnodeprobs = NEW_NUMERIC(m_pProbLists[node]->nProbSize));
			pf = NUMERIC_POINTER(pnodeprobs);
			if (pf && m_pProbLists[node]->pProbs) {
				memcpy(pf, m_pProbLists[node]->pProbs, m_pProbLists[node]->nProbSize * sizeof(double));
			}
			SET_VECTOR_ELT(plist, node, pnodeprobs);
			UNPROTECT(1);
		}
	}
	else { 
		for(node = 0; node < m_numNodes; node++) {
			if (pslotcats) {
				memset(pslotcats, 0, m_maxParents*sizeof(int));
			}
			pnodeprobs = genProbList(node, 0, pslotcats);
			SET_VECTOR_ELT(plist, node, pnodeprobs);
			if(pnodeprobs != R_NilValue)
				UNPROTECT(1);
		}
	}
	SET_SLOT(cnet, install("probabilities"), plist);
	UNPROTECT(1); /* plist */

	if (pslotcats)
		CATNET_FREE(pslotcats);

	pstr = PROTECT(mkString("catNetwork object"));
	SET_SLOT(cnet, install("meta"), pstr);
	UNPROTECT(1); /* pstr */

	PROTECT(pint = NEW_INTEGER(1));
	INTEGER_POINTER(pint)[0] = complexity();
	SET_SLOT(cnet, install("complexity"), pint);
	UNPROTECT(1); /* pint */
	
	PROTECT(pint = NEW_NUMERIC(1));
	floglik = loglik();
	if(floglik > -FLT_MAX)
		NUMERIC_POINTER(pint)[0] = floglik;
	else
		NUMERIC_POINTER(pint)[0] = R_NegInf;

	SET_SLOT(cnet, install("likelihood"), pint);
	UNPROTECT(1); /* pint */

	PROTECT(pint = NEW_INTEGER(m_numNodes));
	INTEGER_POINTER(pint)[0] = m_numNodes;
	for(node = 0; node < m_numNodes; node++) {
		if(m_pProbLists[node])
			INTEGER_POINTER(pint)[node] = m_pProbLists[node]->sampleSize;
	}
	SET_SLOT(cnet, install("nodeSampleSizes"), pint);
	UNPROTECT(1); /* pint */

	UNPROTECT(1); /* cnet */
	return cnet;
}

SEXP RCatnet::genProbList(int node, int paridx, int *pcats) {
	int j, npar;
	SEXP problist;
	double *pslot, *pp;

	if(m_pProbLists == 0 || m_pProbLists[node] == 0 || paridx < 0)
		return R_NilValue;

	if(paridx >= m_numParents[node]) {
		pslot = m_pProbLists[node]->find_slot(0, pcats, 0);
		PROTECT(problist = NEW_NUMERIC(m_numCategories[node]));
		pp = NUMERIC_POINTER(problist);
		if (pp && pslot) {
			memcpy(pp, pslot, m_numCategories[node]*sizeof(double));
		}
		return problist;
	}

	npar = m_parents[node][paridx];
	PROTECT(problist = allocVector(VECSXP, m_numCategories[npar]));
	for(j = 0; j < m_numCategories[npar]; j++) {
		pcats[paridx] = j;
		SET_VECTOR_ELT(problist, j, genProbList(node, paridx + 1, pcats));
		UNPROTECT(1);
	}

	return problist;
}

char *gen_prob_string(int node, SEXP parlist, int paridx, SEXP catlist, SEXP problist, char *str) {
	int j, npar, nlen;
	SEXP parprobs, pcats;
	char *newstr, *aux, *aux2, *aux3;

	if(!str) {
		str = (char*)CATNET_MALLOC(1);
		str[0] = 0;
	}
	if(!str) 
		return 0;

	if(paridx >= length(parlist)) {
		pcats = VECTOR_ELT(catlist, node);
		newstr = (char*)CATNET_MALLOC(((strlen(str)+1+32)*length(pcats))*sizeof(char));
		if (!newstr) {
			CATNET_FREE(str);
			return 0;
		}
		nlen = 0;
		for(j = 0; j < length(pcats); j++) {
			nlen += sprintf(newstr+nlen, "%s%s %f\n", str, CHAR(STRING_ELT(pcats, j)), NUMERIC_POINTER(problist)[j]);
		}
		CATNET_FREE(str);
		str = newstr;
		return str;
	}

	npar = INTEGER_POINTER(parlist)[paridx] - 1;
	pcats = VECTOR_ELT(catlist, npar);

	newstr = (char*)CATNET_MALLOC(sizeof(char));
	if (!newstr) {
		CATNET_FREE(str);
		return 0;
	}

	newstr[0] = 0;
	for(j = 0; j < length(pcats); j++) {
		parprobs = VECTOR_ELT(problist, j);

		aux = (char*)CATNET_MALLOC((strlen(str)+1+8)*sizeof(char));
		if (aux)
			sprintf(aux, "%s%s", str, CHAR(STRING_ELT(pcats, j)));
		aux2 = gen_prob_string(node, parlist, paridx + 1, catlist, parprobs, aux);

		aux3 = (char*)CATNET_MALLOC((strlen(newstr)+strlen(aux2)+2)*sizeof(char));
		if (aux3)
			sprintf(aux3, "%s%s", newstr, aux2);
		CATNET_FREE(newstr);
		newstr = aux3;

		CATNET_FREE(aux2);
	}
	CATNET_FREE(str);
	str = newstr;

	return str;
}

SEXP prob_string(SEXP rnodes, SEXP rparents, SEXP rcatlist, SEXP rproblist) {

	int node;
	SEXP nodepars, nodeproblist, pstr;
	char *str = NULL, *newstr, *aux;

	PROTECT(rnodes = AS_LIST(rnodes));
	PROTECT(rparents = AS_LIST(rparents));
	PROTECT(rcatlist = AS_LIST(rcatlist));
	PROTECT(rproblist = AS_LIST(rproblist));

	for(node = 0; node < length(rnodes); node++) {
		nodepars = VECTOR_ELT(rparents, node);
		nodeproblist = VECTOR_ELT(rproblist, node);
		newstr = gen_prob_string(node, nodepars, 0, rcatlist, nodeproblist, NULL);
		if(str && newstr) {
			aux = (char*)CATNET_MALLOC((strlen(str)+strlen(newstr)+1+16)*sizeof(char));
			if (aux)
				sprintf(aux, "%sNode [%d]:\n%s", str, node, newstr);
			CATNET_FREE(str);
			CATNET_FREE(newstr);
			str = aux;
		}
		else {
			CATNET_FREE(str);
			str = newstr;
		}
	}

	UNPROTECT(4);

	PROTECT(pstr = allocVector(STRSXP, 1));
	SET_STRING_ELT(pstr, 0, mkChar(str));
	UNPROTECT(1);

	return pstr;
}


void gen_prob_vector(int node, SEXP parlist, int paridx, SEXP catlist, SEXP problist, double *&pvec, int &nvec) {
	int j, npar;
	SEXP parprobs, pcats;
	double *newvec;

	if(!pvec) {
		pvec = (double*)CATNET_MALLOC(sizeof(double));
		nvec = 0;
	}
	if(!pvec) 
		return;

	if(paridx >= length(parlist)) {
		pcats = VECTOR_ELT(catlist, node);
		if (length(problist) != length(pcats)) {
			Rprintf("gen_prob_vector: %d:  %d, %d\n", node, length(problist), length(pcats));
			error("Wrong probability table");
			return;
		}
		newvec = (double*)CATNET_MALLOC((nvec + length(pcats))*sizeof(double));
		if (newvec) {
			memcpy(newvec, pvec, nvec*sizeof(double));
			for(j = 0; j < length(pcats); j++) {
				newvec[nvec+j] = NUMERIC_POINTER(problist)[j];
			}
		}
		CATNET_FREE(pvec);
		pvec = newvec;
		nvec += length(pcats);
		return;
	}

	npar = INTEGER_POINTER(parlist)[paridx] - 1;
	pcats = VECTOR_ELT(catlist, npar);
	if (length(problist) != length(pcats)) {
		Rprintf("gen_prob_vector: %d:  %d, %d\n", node, length(problist), length(pcats));
		error("Wrong probability table");
		return;
	}
	for(j = 0; j < length(pcats); j++) {
		parprobs = VECTOR_ELT(problist, j);
		gen_prob_vector(node, parlist, paridx + 1, catlist, parprobs, pvec, nvec);
	}
}

SEXP prob_vector(SEXP rnodes, SEXP rparents, SEXP rcatlist, SEXP rproblist) {

	int node, nvec;
	SEXP nodepars, nodeproblist, pstr;
	SEXP rvec;
	double *pvec, *prvec;

	PROTECT(rnodes = AS_LIST(rnodes));
	PROTECT(rparents = AS_LIST(rparents));
	PROTECT(rcatlist = AS_LIST(rcatlist));
	PROTECT(rproblist = AS_LIST(rproblist));

	PROTECT(pstr = allocVector(VECSXP, length(rnodes)));

	for(node = 0; node < length(rnodes); node++) {
		nodepars = VECTOR_ELT(rparents, node);
		nodeproblist = VECTOR_ELT(rproblist, node);
		pvec = 0;
		nvec = 0;
		gen_prob_vector(node, nodepars, 0, rcatlist, nodeproblist, pvec, nvec);
		PROTECT(rvec = NEW_NUMERIC(nvec));
		prvec = NUMERIC_POINTER(rvec);
		if (prvec && pvec)
			memcpy(prvec, pvec, nvec*sizeof(double));
		CATNET_FREE(pvec);
		SET_VECTOR_ELT(pstr, node, rvec);
		UNPROTECT(1);
	}

	UNPROTECT(5);
	return pstr;
}

SEXP RCatnet::genSamples(SEXP rNumSamples, SEXP rPerturbations, SEXP rNaRate) {

	SEXP rsamples = R_NilValue;
	int numsamples;
	int *pSamples, *pRsamples, *pPerturbations;
	double fNaRate;
	int i, j, k, nnode, *pnodepars, *pnodesample, *porder;
	double u, v, *pnodeprob;
	PROB_LIST<double> *pProbList; 
	
	PROTECT(rNumSamples = AS_INTEGER(rNumSamples));
	numsamples = INTEGER_POINTER(rNumSamples)[0];
	UNPROTECT(1);

	PROTECT(rNaRate = AS_NUMERIC(rNaRate));
	fNaRate = NUMERIC_POINTER(rNaRate)[0];
	UNPROTECT(1);

	porder = getOrder();
	if(!porder)
		return R_NilValue;

	pSamples = (int*)CATNET_MALLOC(m_numNodes*numsamples*sizeof(int));
	if(!pSamples) {
		CATNET_FREE(porder);
		return R_NilValue;
	}

	pnodesample = 0;
	if(m_maxParents > 0) {
		pnodesample = (int*)CATNET_MALLOC(m_maxParents * sizeof(int));
		if(!pnodesample) {
			CATNET_FREE(pSamples);
			CATNET_FREE(porder);
			return R_NilValue;
		}
	}

	pPerturbations = NULL;
	if(!isNull(rPerturbations)) {
		PROTECT(rPerturbations = AS_INTEGER(rPerturbations));
		pPerturbations = INTEGER_POINTER(rPerturbations);
	}

	GetRNGstate();
	for(k = 0; k < m_numNodes; k++) {
		nnode = porder[k];
		pnodepars = m_parents[nnode];
		pProbList = (PROB_LIST<double>*)getNodeProb(nnode);
		if (!pProbList) 
			break;

		for (j = 0; j < numsamples; j++) {

			if(pPerturbations) {
				if(!R_IsNA(pPerturbations[j*m_numNodes + nnode]) && 
					pPerturbations[j*m_numNodes + nnode] >= 1 && pPerturbations[j*m_numNodes + nnode] <= m_numCategories[nnode]) {
					pSamples[j * m_numNodes + nnode] = pPerturbations[j * m_numNodes + nnode];
					continue;
				}
			}

			for (i = 0; i < m_numParents[nnode]; i++) {
				if (pnodepars[i] < 0 || pnodepars[i] >= m_numNodes)
					break;
				pnodesample[i] = (int)(pSamples[j * m_numNodes + pnodepars[i]] - 1);
			}
			pnodeprob = pProbList->find_slot(0, pnodesample, 0);

			u = (double)unif_rand();
			v = 0;
			for(i = 0; i < m_numCategories[nnode]; i++) {
				v += pnodeprob[i];
				if(u <= v)
					break;
			}
			pSamples[j * m_numNodes + nnode] = i + 1;
		}
	}

	k = (int)(fNaRate*m_numNodes);
	if(k > 0 && k < m_numNodes) {
		int ii, fmax, *paux;
		paux = (int*)CATNET_MALLOC(m_numNodes*sizeof(int));
		if (paux) {
			for (j = 0; j < numsamples; j++) {
				for(ii = 0; ii < m_numNodes; ii++)
					paux[ii] = (int)(m_numNodes*unif_rand());
				for(i = 0; i < k; i++) {
					nnode = 0;
					fmax = -(int)RAND_MAX;
					for(ii = 0; ii < m_numNodes; ii++) {
						if(fmax < paux[ii]) {
							fmax = paux[ii];
							nnode = ii;
						}
					}
					if(nnode >= 0 && nnode < m_numNodes) {
						paux[nnode] = -(int)RAND_MAX;
						pSamples[j * m_numNodes + nnode] = R_NaInt;
					}
				}
			}
			CATNET_FREE(paux);
		}
	}
	PutRNGstate();

	if(!isNull(rPerturbations)) {
		UNPROTECT(1); // rPerturbations
	}

	if(pnodesample)
		CATNET_FREE(pnodesample);		
	if(porder)
		CATNET_FREE(porder);

	// output the new matrix
	PROTECT(rsamples = NEW_INTEGER(m_numNodes * numsamples));
	pRsamples = (int*)INTEGER_POINTER(rsamples);
	if (pRsamples && pSamples)
		memcpy(pRsamples, pSamples, m_numNodes*numsamples*sizeof(int));
	UNPROTECT(1);

	CATNET_FREE(pSamples);

	return rsamples;
}
