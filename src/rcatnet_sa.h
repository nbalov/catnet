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
 * rcatnet_sa.h
 *
 *  Created on: Sep 21, 2009
 *      Author: nbalov
 */

/* 
 * version 1.15.1  12dec2016
 */

#ifndef RCATNET_SEARCH_SA_H
#define RCATNET_SEARCH_SA_H

#include "catnet_search2.h"

#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>

#define MAX_NODE_NAME	16

class RCatnetSearchSA {
protected:
	int m_nDrives;
	CATNET_SEARCH2<char, MAX_NODE_NAME, double> **m_pDrives;
	int m_numNodes, m_numSamples, m_maxParentSet;
	int m_bUseCache;

	int *m_pOptOrder;
	int m_nOptNets;
	int **m_pTestOrder, **m_pTestOrderInverse;
	CATNET<char, MAX_NODE_NAME, double> **m_pOptNets;

	void _release();
	int *_genOrder(const int *porder, int norder, int shuffles, int bjump);
	int *_genOrderFormDirProbs(const int *porder, int numnodes, double *matEdgeLiks, double *pOrderProb);

public:
	SEARCH_PARAMETERS **m_pSearchParams;

public:
	RCatnetSearchSA();
	~RCatnetSearchSA();

	SEXP search(SEXP rNodeNames, SEXP rSamples, SEXP rPerturbations, 
			SEXP rMaxParents, SEXP rParentSizes, SEXP rMaxComplexity, SEXP rNodeCats, 
			SEXP rParentsPool, SEXP rFixedParentsPool, SEXP rMaxParentsPool, 
			SEXP rMatEdgeLiks, SEXP rDirProbs, 
			SEXP rModel, SEXP rStartOrder,
			SEXP rTempStart, SEXP rTempCoolFact, SEXP rTempCheckOrders, 
			SEXP rMaxIter, SEXP rOrderShuffles, SEXP rStopDiff, 
			SEXP rThreads, SEXP rUseCache, SEXP rEcho);

};

#endif /* RCATNET_SEARCH_SA_H */
