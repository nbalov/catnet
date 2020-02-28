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
 * catnet_rexport.cpp
 *
 *  Created on: Sep 21, 2009
 *      Author: nbalov
 */

/* 
 * version 1.15.6  25feb2020
 */

#include "utils.h"
#include "rcatnet.h"
#include "rcatnet_search.h"
#include "rcatnet_sa.h"
#include "rcatnet_hist.h"

extern "C" {

extern int g_setseed;
extern size_t g_memcounter;

SEXP catnetReleaseCache()
{
	ReleaseCache();
	return R_NilValue;
}

SEXP createRCatnet(SEXP cnet)
{
	SEXP pcnet = R_NilValue;
	PROTECT(cnet);
	RCatnet *rnet = new RCatnet(cnet);
	UNPROTECT(1);
	pcnet = rnet->genRcatnet((const char*)"catNetwork");
	delete rnet;
	return pcnet;
}

SEXP catnetMarginalProb(SEXP cnet, SEXP rnode)
{
	int node, i, ncats;
	SEXP rvec = R_NilValue;
	double *pvec;

	if(!isInteger(AS_INTEGER(rnode)))
		error("node should be an integer");
	PROTECT(rnode = AS_INTEGER(rnode));
	node = INTEGER_POINTER(rnode)[0];
	UNPROTECT(1);

	PROTECT(cnet);
	RCatnet *rnet = new RCatnet(cnet);
	UNPROTECT(1);

	if(!rnet)
		return rvec;

	if(node < 1 || node > rnet->numNodes())
		return rvec;
	node--;

	double *pmarg = rnet->marginal_prob(node);
	if(!pmarg)
		return rvec;

	ncats = rnet->numCategories(node);
	PROTECT(rvec = NEW_NUMERIC(ncats));
	pvec = NUMERIC_POINTER(rvec);
	for(i = 0; i < ncats; i++) {
		pvec[i] = pmarg[i];
	}
	UNPROTECT(1);

	CATNET_FREE(pmarg);
	delete rnet;

	//char str[128];
	//sprintf(str, "Mem Balance  %d\n", (int)g_memcounter);
	//printf(str);

	return rvec;
}

SEXP catnetJointProb(SEXP cnet, SEXP rnode)
{
	int node, jointprobsize;
	SEXP rvec = R_NilValue;
	double *pvec;

	PROTECT(rnode = AS_INTEGER(rnode));
	node = INTEGER_VALUE(rnode);
	UNPROTECT(1);

	PROTECT(cnet);
	RCatnet *rnet = new RCatnet(cnet);
	UNPROTECT(1);

	if(!rnet)
		return rvec;
	if(node < 1 || node > rnet->numNodes())
		return rvec;
	node--;

	jointprobsize = 0;
	double *pjoint = rnet->findJointProb(node, jointprobsize);
	if(!pjoint) {
		delete rnet;
		return rvec;
	}

	PROTECT(rvec = NEW_NUMERIC(jointprobsize));
	pvec = NUMERIC_POINTER(rvec);
	if (pvec && pjoint)
		memcpy(pvec, pjoint, jointprobsize*sizeof(double));
	UNPROTECT(1);

	CATNET_FREE(pjoint);

	delete rnet;
	return rvec;
}

SEXP catnetFindParentPool(SEXP cnet, SEXP rnode)
{
	int node, i, poolsize;
	SEXP rvec = R_NilValue;
	int *ppool, *pvec;

	PROTECT(rnode = AS_INTEGER(rnode));
	node = INTEGER_VALUE(rnode);
	UNPROTECT(1);

	PROTECT(cnet);
	RCatnet *rnet = new RCatnet(cnet);
	UNPROTECT(1);

	if(!rnet)
		return rvec;

	if(node < 1 || node > rnet->numNodes())
		return rvec;
	node--;

	ppool = rnet->findParentPool(node, poolsize);
	if (!ppool) {
		delete rnet;
		return rvec;
	}

	PROTECT(rvec = NEW_INTEGER(poolsize));
	pvec = INTEGER_POINTER(rvec);
	for(i = 0; i < poolsize; i++) {
		pvec[i] = ppool[i]+1;
	}
	UNPROTECT(1);

	CATNET_FREE(ppool);

	delete rnet;
	return rvec;
}

SEXP show_catnet(SEXP rnodes, SEXP rparents, SEXP rcatlist, SEXP rproblist)
{
	int i, nlen, m_numNodes, nnode;
	char *strbuff;
	SEXP pf, pstr;

	m_numNodes = length(rnodes);
	strbuff = (char*)CATNET_MALLOC(16+m_numNodes*m_numNodes*(2+MAX_NODE_NAME));
	if (!strbuff) {
		return R_NilValue;
	}

	PROTECT(rnodes = AS_LIST(rnodes));
	PROTECT(rparents = AS_LIST(rparents));
	PROTECT(rcatlist = AS_LIST(rcatlist));
	PROTECT(rproblist = AS_LIST(rproblist));
	PROTECT(pstr = allocVector(STRSXP, 3));

	nlen = sprintf(strbuff, "Nodes = %d: ", m_numNodes);
	for(nnode = 0; nnode < m_numNodes; nnode++) {
		PROTECT(pf = VECTOR_ELT(rnodes, nnode));
		if(IS_VECTOR(pf)) {
			nlen += sprintf(strbuff+nlen, "%s, ", CHAR(STRING_ELT(pf, 0)));
		}
		UNPROTECT(1);
	}
	nlen += sprintf(strbuff+nlen, "\n");
	SET_STRING_ELT(pstr, 0, mkChar(strbuff));

	nlen = sprintf(strbuff, "Parents:\n");
	for(nnode = 0; nnode < m_numNodes; nnode++) {
		PROTECT(pf = VECTOR_ELT(rparents, nnode));
		nlen += sprintf(strbuff+nlen, "[%d] ", nnode);
		if(IS_VECTOR(pf)) {
			for(i = 0; i < length(pf); i++)
				nlen += sprintf(strbuff+nlen, "%d, ", INTEGER_POINTER(pf)[i]-1);
			nlen += sprintf(strbuff+nlen, "\n");
		}
		else 
			nlen += sprintf(strbuff+nlen, "\n");
		UNPROTECT(1);
	}
	SET_STRING_ELT(pstr, 1, mkChar(strbuff));

	nlen = sprintf(strbuff, "Categories:\n");
	for(nnode = 0; nnode < m_numNodes; nnode++) {
		PROTECT(pf = VECTOR_ELT(rcatlist, nnode));
		if(IS_VECTOR(pf)) {
			for(i = 0; i < length(pf); i++)
				nlen += sprintf(strbuff+nlen, "%s, ", CHAR(STRING_ELT(pf,i)));
			nlen += sprintf(strbuff+nlen, "\n");
		}
		UNPROTECT(1);
	}
	SET_STRING_ELT(pstr, 2, mkChar(strbuff));

	UNPROTECT(5);

	CATNET_FREE(strbuff);

	return pstr;
}

SEXP showCatnet(SEXP cnet)
{
	SEXP rnodes, rparents, rcatlist, rproblist;

	PROTECT(cnet);
	//if(isS4(cnet))
	//	printf("cnet is an object\n");

	rnodes = GET_SLOT(cnet, install("nodes"));
	//if(IS_VECTOR(m_nodeNames))
	//	printf("m_nodeNames is a vector\n");
	rparents = GET_SLOT(cnet, install("parents"));
	rcatlist = GET_SLOT(cnet, install("categories"));
	rproblist = GET_SLOT(cnet, install("probabilities"));

	if(rnodes == R_NilValue || rparents == R_NilValue || rcatlist == R_NilValue || rproblist == R_NilValue) {
		UNPROTECT(1);
		return R_NilValue;
	}

	SEXP res = show_catnet(rnodes, rparents, rcatlist, rproblist);

	UNPROTECT(1);

	return res;
}

SEXP catnetOptimalNetsForOrder(SEXP rSamples, SEXP rPerturbations, 
	SEXP rMaxParents, SEXP rParentSizes, SEXP rMaxComplexity, 
	SEXP rOrder, SEXP rNodeCats, 
	SEXP rParentsPool, SEXP rFixedParentsPool, SEXP rMatEdgeLiks, 
	SEXP rUseCache, SEXP rEcho) {

	//if(!isMatrix(rSamples))
	//	error("Data is not a matrix");
	if(!isNull(rPerturbations) && !isMatrix(rPerturbations))
		error("Perturbations should be a matrix");
	if(!isInteger(AS_INTEGER(rMaxParents)))
		error("maxParents should be an integer");
	if(!isNull(rParentSizes) && !isVector(rParentSizes))
		error("ParentSizes should be a vector");
	if(!isInteger(AS_INTEGER(rMaxComplexity)))
		error("maxComplexity should be an integer");
	if(!isVector(rOrder))
		error("Order should be a vector");
	if(!isNull(rNodeCats) && !isVector(rNodeCats))
		error("NodeCats should be a list");
	if(!isNull(rParentsPool) && !isVector(rParentsPool))
		error("ParentsPool should be a list");
	if(!isNull(rFixedParentsPool) && !isVector(rFixedParentsPool))
		error("FixedParentsPool should be a list");
	if(!isNull(rMatEdgeLiks) && !isMatrix(rMatEdgeLiks))
		error("rMatEdgeLiks should be a matrix");
	if(!isNull(rUseCache) && !isLogical(rUseCache))
		error("UseCache should be logical");
	if(!isNull(rEcho) && !isLogical(rEcho))
		error("Echo should be logical");

	RCatnetSearch * pengine = new RCatnetSearch;
	SEXP res = pengine -> estimateCatnets(rSamples, rPerturbations, 
						rMaxParents, rParentSizes, 
						rMaxComplexity, rOrder, rNodeCats, 
				                rParentsPool, rFixedParentsPool, rMatEdgeLiks, 
						rUseCache, rEcho);
	delete pengine;

	//char str[128];
	//sprintf(str, "Mem Balance  %d\n", (int)g_memcounter);
	//printf(str);

	return res;
}

SEXP catnetOptimalNetsSA(SEXP rNodeNames, SEXP rSamples, SEXP rPerturbations, 
	SEXP rMaxParents, SEXP rParentSizes, SEXP rMaxComplexity, SEXP rNodeCats, 
	SEXP rParentsPool, SEXP rFixedParentsPool, SEXP rMaxParentsPool, 
	SEXP rMatEdgeLiks, SEXP rDirProbs, 
	SEXP rModel, SEXP rStartOrder,
	SEXP rTempStart, SEXP rTempCoolFact, SEXP rTempCheckOrders, 
	SEXP rMaxIter, SEXP rOrderShuffles, SEXP rStopDiff, 
	SEXP rThreads, SEXP rUseCache, SEXP rEcho) {

	//if(!isMatrix(rSamples))
	//	error("Data is not a matrix");
	if(!isNull(rPerturbations) && !isMatrix(rPerturbations))
		error("Perturbations should be a matrix");
	if(!isInteger(AS_INTEGER(rMaxParents)))
		error("maxParents should be an integer");
	if(!isNull(rParentSizes) && !isVector(rParentSizes))
		error("ParentSizes should be a vector");
	if(!isInteger(AS_INTEGER(rMaxComplexity)))
		error("maxComplexity should be an integer");
	if(!isVector(rStartOrder))
		error("startOrder should be a vector");
	if(!isNull(rNodeCats) && !isVector(rNodeCats))
		error("NodeCats should be a list");
	if(!isNull(rParentsPool) && !isVector(rParentsPool))
		error("ParentsPool should be a list");
	if(!isNull(rFixedParentsPool) && !isVector(rFixedParentsPool))
		error("FixedParentsPool should be a list");
	if(!isInteger(AS_INTEGER(rMaxParentsPool)))
		error("maxParentsPool should be an integer");
	if(!isNull(rMatEdgeLiks) && !isMatrix(rMatEdgeLiks))
		error("rMatEdgeLiks should be a matrix");
	if(!isNull(rDirProbs) && !isMatrix(rDirProbs))
		error("rDirProbs should be a matrix");
	if(!isNumeric(AS_NUMERIC(rTempStart)))
		error("tempStart should be numerical");
	if(!isNumeric(AS_NUMERIC(rTempCoolFact)))
		error("coolFact should be numerical");
	if(!isNumeric(AS_NUMERIC(rTempCheckOrders)))
		error("tempCheckOrders should be numerical");
	if(!isInteger(AS_INTEGER(rMaxIter)))
		error("maxIter should be an integer");
	if(!isNumeric(AS_NUMERIC(rStopDiff)))
		error("stopDiff should be numerical");
	if(!isNumeric(AS_NUMERIC(rOrderShuffles)))
		error("orderShuffles should be numerical");
	if(!isInteger(AS_INTEGER(rThreads)))
		error("Threads should be an integer");
	if(!isNull(rUseCache) && !isLogical(rUseCache))
		error("UseCache should be logical");
	if(!isNull(rEcho) && !isLogical(rEcho))
		error("Echo should be logical");

	RCatnetSearchSA * pengine = new RCatnetSearchSA;
	SEXP res = pengine -> search(rNodeNames, rSamples, rPerturbations, 
			rMaxParents, rParentSizes, rMaxComplexity, rNodeCats, 
			rParentsPool, rFixedParentsPool, rMaxParentsPool, 
			rMatEdgeLiks, rDirProbs, 
			rModel, rStartOrder,
			rTempStart, rTempCoolFact, rTempCheckOrders, 
			rMaxIter, rOrderShuffles, rStopDiff, 
			rThreads, rUseCache, rEcho);
	delete pengine;

	//char str[128];
	//sprintf(str, "Mem Balance %d\n", (int)g_memcounter);
	//printf(str);

	return res;
}

SEXP catnetParHistogram(SEXP rSamples, SEXP rPerturbations, 
	SEXP rMaxParents, SEXP rParentSizes, 
	SEXP rMaxComplexity, SEXP rNodeCats, 
	SEXP rParentsPool, SEXP rFixedParentsPool, 
	SEXP rScore, SEXP rWeight, SEXP rMaxIter,
	SEXP rThreads, SEXP rUseCache, SEXP rEcho)
{
	//if(!isMatrix(rSamples))
	//	error("Data is not a matrix");
	if(!isNull(rPerturbations) && !isMatrix(rPerturbations))
		error("Perturbations should be a matrix");
	if(!isInteger(AS_INTEGER(rMaxParents)))
		error("maxParents should be an integer");
	if(!isNull(rParentSizes) && !isVector(rParentSizes))
		error("ParentSizes should be a vector");
	if(!isInteger(AS_INTEGER(rMaxComplexity)))
		error("maxComplexity should be an integer");
	if(!isNull(rNodeCats) && !isVector(rNodeCats))
		error("NodeCats should be a list");
	if(!isNull(rParentsPool) && !isVector(rParentsPool))
		error("ParentsPool should be a list");
	if(!isNull(rFixedParentsPool) && !isVector(rFixedParentsPool))
		error("FixedParentsPool should be a list");
	if(!isInteger(AS_INTEGER(rMaxIter)))
		error("maxIter should be an integer");
	if(!isInteger(AS_INTEGER(rWeight)))
		error("weight should be an integer");
	if(!isInteger(AS_INTEGER(rThreads)))
		error("Threads should be an integer");
	if(!isNull(rUseCache) && !isLogical(rUseCache))
		error("UseCache should be logical");
	if(!isNull(rEcho) && !isLogical(rEcho))
		error("Echo should be logical");

	RCatnetSearchHist * pengine = new RCatnetSearchHist;
	SEXP res = pengine -> search(rSamples, rPerturbations, 
			rMaxParents, rParentSizes, rMaxComplexity, rNodeCats, 
			rParentsPool, rFixedParentsPool, 
			rScore, rWeight, rMaxIter, 
			rThreads, rUseCache, rEcho);
	delete pengine;

	//char str[128];
	//sprintf(str, "Mem Balance  %d\n", (int)g_memcounter);
	//printf(str);

	return res;
}

SEXP catnetSetProb(SEXP cnet, SEXP rSamples, SEXP rPerturbations) {

	int *pSamples, *psubSamples, *pPerturbations;
	int nnode, numsamples, numnodes, numsubsamples, j;
	SEXP dim;
	SEXP pcnet = R_NilValue;

	if(!isMatrix(rSamples))
		error("Data is not a matrix");
	if(!isNull(rPerturbations) && !isMatrix(rPerturbations))
		error("Perturbations is not a vector");

	PROTECT(cnet);
	RCatnet *rnet = new RCatnet(cnet);
	UNPROTECT(1);

	PROTECT(rSamples = AS_INTEGER(rSamples));
	pSamples = INTEGER(rSamples);

	dim = GET_DIM(rSamples);
	numnodes = INTEGER(dim)[0];
	numsamples = INTEGER(dim)[1];

	////////////////////////////////////////
	// Danger Ahead
	// We don's check that sample nodes actually correspond to the cnet's nodes
	// Missmatch of categories possible

	// pSamples are assumed positive indices
	for(j = 0; j < numnodes*numsamples; j++) {
		if(R_IsNA(pSamples[j]) || pSamples[j] < 1)
			pSamples[j] = CATNET_NAN;
		else
			pSamples[j]--;
	}

	psubSamples = 0;
	pPerturbations = 0;
	if(!isNull(rPerturbations)) {
		PROTECT(rPerturbations = AS_INTEGER(rPerturbations));
		pPerturbations = INTEGER_POINTER(rPerturbations);
		psubSamples = (int*)CATNET_MALLOC(numnodes*numsamples*sizeof(int));
	}

	for(nnode = 0; nnode < numnodes; nnode++) {
		if(pPerturbations && psubSamples) {
			numsubsamples = 0;
			for(j = 0; j < numsamples; j++) {		
				if(!pPerturbations[j * numnodes + nnode]) {
					memcpy(psubSamples + numsubsamples*numnodes, 
						pSamples + j*numnodes, numnodes*sizeof(int));
					numsubsamples++;
				}
			}
			rnet->setNodeSampleProb(nnode, psubSamples, numsubsamples, 1);
		}
		else {
			rnet->setNodeSampleProb(nnode, pSamples, numsamples, 1);
		}
	}

	UNPROTECT(1);
	if(pPerturbations) {
		UNPROTECT(1);
	}

	if(psubSamples)
		CATNET_FREE(psubSamples);
	pcnet = rnet->genRcatnet((const char*)"catNetwork");
	delete rnet;

	//char str[128];
	//sprintf(str, "Mem Balance  %d\n", (int)g_memcounter);
	//printf(str);

	return pcnet;

}

SEXP catnetLoglik(SEXP cnet, SEXP rSamples, SEXP rPerturbations, SEXP rBySample) {

	int *pSamples, *pPerturbations;
	int numsamples, numnodes, j, bysample;
	double *floglik, *pvec;
	SEXP dim, rvec = R_NilValue;

	if(!isMatrix(rSamples))
		error("Data should be a matrix");
	if(!isNull(rPerturbations) && !isMatrix(rPerturbations))
		error("Perturbations should be a matrix");

	bysample = 0;
	PROTECT(rBySample = AS_LOGICAL(rBySample));
	bysample = LOGICAL(rBySample)[0];
	UNPROTECT(1);

	////////////////////////////////////////
	// Danger Ahead
	// We don's check that sample nodes actually correspond to the cnet's nodes
	// Missmatch of categories possible

	PROTECT(cnet);
	RCatnet *rnet = new RCatnet(cnet);
	UNPROTECT(1);

	PROTECT(rSamples = AS_INTEGER(rSamples));
	pSamples = INTEGER(rSamples);

	dim = GET_DIM(rSamples);
	numnodes = INTEGER(dim)[0];
	numsamples = INTEGER(dim)[1];

	// pSamples are assumed positive indices
	for(j = 0; j < numnodes*numsamples; j++) {
		if(R_IsNA(pSamples[j]) || pSamples[j] < 1)
			pSamples[j] = CATNET_NAN;
		else
			pSamples[j]--;
	}

	pPerturbations = 0;
	if(!isNull(rPerturbations)) {
		PROTECT(rPerturbations = AS_INTEGER(rPerturbations));
		pPerturbations = INTEGER_POINTER(rPerturbations);
	}

	if(bysample)
		floglik = rnet->bySampleLoglikVector(pSamples, numsamples, pPerturbations);
	else
		floglik = rnet->sampleLoglikVector(pSamples, numsamples, pPerturbations);

	UNPROTECT(1);
	delete rnet;

	if(pPerturbations)
		UNPROTECT(1);

	if(floglik) {
		if(bysample) {
			PROTECT(rvec = NEW_NUMERIC(numsamples));
			pvec = NUMERIC_POINTER(rvec);
			for(j = 0; j < numsamples; j++) {
				pvec[j] =  R_NegInf;
				if(floglik[j] > -FLT_MAX)
					pvec[j] = floglik[j];
			}
		}
		else {
			PROTECT(rvec = NEW_NUMERIC(numnodes));
			pvec = NUMERIC_POINTER(rvec);
			for(j = 0; j < numnodes; j++) {
				pvec[j] =  R_NegInf;
				if(floglik[j] > -FLT_MAX)
					pvec[j] = floglik[j];
			}
		}
		UNPROTECT(1);
		CATNET_FREE(floglik);
	}

	//char str[128];
	//sprintf(str, "Mem Balance  %d\n", (int)g_memcounter);
	//printf(str);

	return rvec;

}


SEXP catnetNodeLoglik(SEXP cnet, SEXP rNode, SEXP rSamples, SEXP rPerturbations) {

	int *pSamples, *pPerturbations;
	int *psubSamples, numsubsamples;
	int numsamples, numnodes, i, j, nnode, nnodes, *pnodes;
	double floglik, *pvec;
	SEXP dim, rvec = R_NilValue;

	if(!isMatrix(rSamples))
		error("Data should be a matrix");
	if(!isNull(rPerturbations) && !isMatrix(rPerturbations))
		error("Perturbations should be a matrix");
	if(!isInteger(AS_INTEGER(rNode)))
		error("Node should be an integer");

	////////////////////////////////////////
	// Danger Ahead
	// We don's check that sample nodes actually correspond to the cnet's nodes
	// Missmatch of categories possible
	PROTECT(cnet);
	RCatnet *rnet = new RCatnet(cnet);
	UNPROTECT(1);

	nnodes = length(rNode);
	if(nnodes < 1)
		return rvec;
	pnodes = (int*)CATNET_MALLOC(nnodes*sizeof(int));
	PROTECT(rNode = AS_INTEGER(rNode));
	if (pnodes && INTEGER_POINTER(rNode))
		memcpy(pnodes, INTEGER_POINTER(rNode), nnodes*sizeof(int));
	UNPROTECT(1);

	PROTECT(rSamples = AS_INTEGER(rSamples));
	pSamples = INTEGER(rSamples);

	dim = GET_DIM(rSamples);
	numnodes = INTEGER(dim)[0];
	numsamples = INTEGER(dim)[1];

	// pSamples are assumed positive indices
	for(j = 0; j < numnodes*numsamples; j++) {
		if(R_IsNA(pSamples[j]) || pSamples[j] < 1)
			pSamples[j] = CATNET_NAN;
		else
			pSamples[j]--;
	}
	PROTECT(rvec = NEW_NUMERIC(nnodes));
	pvec = NUMERIC_POINTER(rvec);

	for(i = 0; i < nnodes; i++) { 
		nnode = pnodes[i] - 1;
		psubSamples = 0;
		pPerturbations = 0;
		floglik = -FLT_MAX;
		if(!isNull(rPerturbations)) {
			PROTECT(rPerturbations = AS_INTEGER(rPerturbations));
			pPerturbations = INTEGER_POINTER(rPerturbations);
			psubSamples = (int*)CATNET_MALLOC(numnodes*numsamples*sizeof(int));
			if (psubSamples) {
				numsubsamples = 0;
				for(j = 0; j < numsamples; j++) {
					if(!pPerturbations[j * numnodes + nnode]) {
						memcpy(psubSamples + numsubsamples*numnodes, pSamples + j*numnodes, numnodes*sizeof(int));
						numsubsamples++;
					}
				}
				floglik = rnet->sampleNodeLoglik(nnode, psubSamples, numsubsamples);
				UNPROTECT(1);
				CATNET_FREE(psubSamples);
			}
		}
		else
			floglik = rnet->sampleNodeLoglik(nnode, pSamples, numsamples);

		pvec[i] = R_NegInf;
		if(floglik > -FLT_MAX)
			pvec[i] =  floglik;
	}
	UNPROTECT(1); // rSamples
	UNPROTECT(1); // rvec

	delete rnet;
	CATNET_FREE(pnodes);

	//char str[128];
	//sprintf(str, "Mem Balance  %d\n", (int)g_memcounter);
	//Rprintf(str);

	return rvec;

}

SEXP catnetSamples(SEXP cnet, SEXP rNumSamples, SEXP rPerturbations, SEXP rNaRate) {

	SEXP rsamples = R_NilValue;

	if(!isInteger(rNumSamples))
		error("rNumSamples should be integer");
	if(!isNumeric(AS_NUMERIC(rNaRate)))
		error("rNaRate should be numerical between 0 and 1");

	PROTECT(cnet);
	RCatnet *rnet = new RCatnet(cnet);
	UNPROTECT(1);
	if(!rnet)
		return R_NilValue;

	rsamples = rnet->genSamples(rNumSamples, rPerturbations, rNaRate);

	delete rnet;

	return rsamples;
}


SEXP catnetSetSeed(SEXP rSeed)
{
	int nSeed;

	if(!isInteger(AS_INTEGER(rSeed)))
		error("The seed should be an integer");

	PROTECT(rSeed = AS_INTEGER(rSeed));
	nSeed = INTEGER_POINTER(rSeed)[0];
	UNPROTECT(1);

	g_setseed = nSeed;

	return(R_NilValue);
}

} // extern "C"
