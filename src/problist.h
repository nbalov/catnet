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
 * problist.h
 *
 *  Created on: Nov 16, 2009
 *      Author: Nikolay Balov
 */

/* 
 * version 1.15.1  12dec2016
 */

#ifndef PROBLIST_H_
#define PROBLIST_H_

using namespace std;

template<class t_prob>
struct PROB_LIST {
	t_prob *pProbs;
	int nProbSize;
	int numCats;
	int numPars;
	int *numParCats;
	int *pBlockSize;
	t_prob loglik;
	t_prob priorlik;
	int sampleSize;

	void reset() {
		if (numParCats)
			CATNET_FREE(numParCats);
		if (pBlockSize)
			CATNET_FREE(pBlockSize);
		if (pProbs)
			CATNET_FREE(pProbs);
		numPars = 0;
		numCats = 0;
		numParCats = 0;
		pBlockSize = 0;
		pProbs = 0;
		nProbSize = 0;
		loglik = 0;
		priorlik = 0;
		sampleSize = 0;
	}


	PROB_LIST() {
		numPars = 0;
		numCats = 0;
		numParCats = 0;
		pBlockSize = 0;
		pProbs = 0;
		nProbSize = 0;
		loglik = 0;
		priorlik = 0;
		sampleSize = 0;
	}

	PROB_LIST<t_prob>& operator =(const PROB_LIST<t_prob> &plist) {
		numPars = plist.numPars;
		numCats = plist.numCats;
		if (numParCats)
			CATNET_FREE(numParCats);
		numParCats = 0;
		if (pBlockSize)
			CATNET_FREE(pBlockSize);
		pBlockSize = 0;
		if(numPars > 0) {
			numParCats = (int*) CATNET_MALLOC(numPars * sizeof(int));
			if (numParCats) {
				memset(numParCats, 0, numPars * sizeof(int));
				if (plist.numParCats)
					memcpy(numParCats, plist.numParCats, numPars * sizeof(int));
			}
			pBlockSize = (int*) CATNET_MALLOC(numPars * sizeof(int));
			if (pBlockSize) {
				memset(pBlockSize, 0, numPars * sizeof(int));
				if (plist.pBlockSize)
					memcpy(pBlockSize, plist.pBlockSize, numPars * sizeof(int));
			}
		}
		nProbSize = plist.nProbSize;
		if (pProbs)
			CATNET_FREE(pProbs);
		pProbs = (t_prob*) CATNET_MALLOC(nProbSize * sizeof(t_prob));
		memset(pProbs, 0, nProbSize * sizeof(t_prob));
		if (pProbs) {
			if (plist.pProbs && nProbSize > 0) {
				for(int i = 0; i < nProbSize; i++)
					pProbs[i] = plist.pProbs[i];
			}
		}
		loglik = plist.loglik;
		priorlik = plist.priorlik;
		sampleSize = plist.sampleSize;
		return *this;
	}

	PROB_LIST(int ncats, int nmaxcats = 2, int npars = 0, int *parcats = 0,
			t_prob *pprobs = 0, int probsize = 0, int samples = 0) {
		int i;
		if (ncats < 1 || nmaxcats < 1 || npars < 0 || (npars && !parcats))
			return;
		numPars = npars;
		numCats = ncats;
		numParCats = 0;
		pBlockSize = 0;
		pProbs = 0;
		nProbSize = 0;
		loglik = 0;
		priorlik = 0;
		sampleSize = samples;
		if (numPars > 0) {
			numParCats = (int*) CATNET_MALLOC(numPars * sizeof(int));
			if (numParCats) {
				if (parcats)
					memcpy(numParCats, parcats, numPars * sizeof(int));
				else
					for (i = 0; i < numPars; i++)
						numParCats[i] = nmaxcats;
			}
			pBlockSize = (int*) CATNET_MALLOC(numPars * sizeof(int));
			if (pBlockSize) {
				pBlockSize[numPars - 1] = ncats;
				for (i = numPars - 1; i > 0; i--) {
					if (parcats[i] < 1 || parcats[i] > nmaxcats) {
						CATNET_FREE(pBlockSize);
						pBlockSize = 0;
						numPars = 0;
						return;
					}
					pBlockSize[i - 1] = parcats[i] * pBlockSize[i];
				}
				nProbSize = pBlockSize[0] * parcats[0];
			}
		}
		else {
			nProbSize = ncats;
		}

		pProbs = (t_prob*) CATNET_MALLOC(nProbSize * sizeof(t_prob));
		if (pProbs) {
			memset(pProbs, 0, nProbSize * sizeof(t_prob));
			if (pProbs && pprobs) {
				if (probsize != nProbSize) {
					return;
				}
				memcpy(pProbs, pprobs, nProbSize * sizeof(t_prob));
			}
		}
	}

	~PROB_LIST() {
		reset();
	}

	void set_zero() {
		if (pProbs)
			memset(pProbs, 0, nProbSize * sizeof(t_prob));
		loglik = 0;
		priorlik = 0;
	}

	void printf() {
		int i;
		Rprintf("pProbs = ");
		for (i = 0; i < nProbSize; i++) {
			Rprintf("%f  ", pProbs[i]);
		}
		Rprintf("\n");
	}

	void normalize() {
		int i, k;
		t_prob pp;
		if (!pProbs)
			return;
		k = 0;
		while (k < nProbSize) {
			pp = 0;
			for (i = 0; i < numCats; i++) {
				pp += pProbs[k + i];
			}
			if (pp > 0) {
				pp = 1 / pp;
				for (i = 0; i < numCats; i++)
					pProbs[k + i] *= pp;
			}
			else {
				// set uniform probability
				t_prob favg;
				favg = (t_prob)1 / (t_prob)numCats;
				pp = 0;
				for (i = 0; i < numCats - 1; i++) {
					pProbs[k + i] = favg;
					pp += favg;			
				}
				pProbs[k + numCats - 1] = 1 - pp;
			}
			k += numCats;
		}
	}

	t_prob loglikelihood() {
		int i, k;
		t_prob pp;
		loglik = 0;
		k = 0;
		while (k < nProbSize) {
			pp = 0;
			for (i = 0; i < numCats; i++)
				pp += pProbs[k + i];
			if (pp <= 0) {
				/////////////////////////////////////////////////////////////////
				// The next condition determines whether parent sets with  
				// non-sample-populated slots should be skipped
				// loglik = -FLT_MAX;
				// break;
			}
			else {
				pp = 1 / pp;
				for (i = 0; i < numCats; i++) {
					if(pProbs[k + i] > 0)
						loglik += pProbs[k + i] * log(pProbs[k + i] * pp);
				}
			}
			k += numCats;
		}
		return(loglik + priorlik); 
	}

	t_prob *find_slot(t_prob* prob, int *pcats, int parid) {
		if (prob == 0)
			prob = pProbs;
		if (parid >= numPars || pcats == 0)
			return prob;
		int parcat = (int)pcats[parid];
		if (parcat < 0 || parcat >= numParCats[parid])
			return 0;
		if (parid == numPars - 1)
			return prob + parcat * pBlockSize[parid];
		return find_slot(prob + parcat * pBlockSize[parid], pcats, parid + 1);
	}

};

#endif /* PROBLIST_H_ */

