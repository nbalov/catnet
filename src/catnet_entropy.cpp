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
 * rcatnet_entropy.cpp
 *
 *  Created on: Sep 21, 2009
 *      Author: nbalov
 */

/* 
 * version 1.15.6  25feb2020
 */

#include "utils.h"
#include "rcatnet.h"

extern "C" {

extern size_t g_memcounter;

SEXP catnetEntropyPairwise(SEXP rSamples, SEXP rPerturbations) {

	int *pSamples, *pPerturbations;
	int *psubSamples, numsubsamples, *psamples1, numsamples1;
	int *pNodeNumCats, **pNodeCats, mincat, maxcat, maxCategories, *pprobs;
	int numsamples, numnodes, i, j, k, d, nnode1, nnode2;
	double floglik, faux, fsum, *pvec, *klmat;
	SEXP dim, rvec = R_NilValue;

	if(!isMatrix(rSamples))
		error("Data should be a matrix");
	if(!isNull(rPerturbations) && !isMatrix(rPerturbations))
		error("Perturbations should be a matrix");

	if(isNull(rSamples)) {
		return rvec;
	}
	PROTECT(rSamples = AS_INTEGER(rSamples));
	pSamples = INTEGER(rSamples);
	
	dim = GET_DIM(rSamples);
	numnodes = INTEGER(dim)[0];
	numsamples = INTEGER(dim)[1];

	// pSamples are assumed positive indices
	for(j = 0; j < numnodes*numsamples; j++) {
		pSamples[j]--;
	}

	// categoies
	pNodeNumCats = (int*)CATNET_MALLOC(numnodes*sizeof(int));
	if (!pNodeNumCats) {
		UNPROTECT(1); //rSamples
		return rvec;
	}
	pNodeCats = (int**)CATNET_MALLOC(numnodes*sizeof(int*));
	if (!pNodeCats) {
		CATNET_FREE(pNodeNumCats);
		UNPROTECT(1); //rSamples
		return rvec;
	}
	memset(pNodeCats, 0, numnodes*sizeof(int*));
	memset(pNodeNumCats, 0, numnodes*sizeof(int));

	maxCategories = 1;
	for(i = 0; i < numnodes; i++) {
		mincat = INT_MAX;
		maxcat = -INT_MAX;
		for(j = 0; j < numsamples; j++) {
			if(pSamples[j*numnodes + i] < mincat)
				mincat = pSamples[j*numnodes + i];
			if(pSamples[j*numnodes + i] > maxcat)
				maxcat = pSamples[j*numnodes + i];
		}
		pNodeNumCats[i] = maxcat - mincat + 1;
		pNodeCats[i] = (int*)CATNET_MALLOC(pNodeNumCats[i]*sizeof(int));
		for(j = 0; j < pNodeNumCats[i]; j++)
			pNodeCats[i][j] = mincat + j;
	}
	for(i = 0; i < numnodes; i++) {
		/* order pNodeNumCats[i] */
		for(j = 0; j < pNodeNumCats[i]; j++) {
			for(k = j + 1; k < pNodeNumCats[i]; k++) {
				if(pNodeCats[i][j] > pNodeCats[i][k]) {
					d = pNodeCats[i][j]; 
					pNodeCats[i][j] = pNodeCats[i][k];
					pNodeCats[i][k] = d;
				}
			}
		} 
		for(j = 0; j < numsamples; j++) {
			for(d = 0; d < pNodeNumCats[i]; d++)
				if(pNodeCats[i][d] == pSamples[j*numnodes + i])
					break;
			pSamples[j*numnodes + i] = d;
		}
		if(maxCategories < pNodeNumCats[i])
			maxCategories = pNodeNumCats[i];
	}

	pprobs = (int*)CATNET_MALLOC(maxCategories*maxCategories*sizeof(int));
	if (!pprobs) {
		CATNET_FREE(pNodeNumCats);
		CATNET_FREE(pNodeCats);
		UNPROTECT(1); //rSamples
		return rvec;
	}

	klmat = (double*)CATNET_MALLOC(numnodes*numnodes*sizeof(double));
	if (!klmat) {
		CATNET_FREE(pNodeNumCats);
		CATNET_FREE(pNodeCats);
		CATNET_FREE(pprobs);
		UNPROTECT(1); //rSamples
		return rvec;
	}
	
	memset(klmat, 0, numnodes*numnodes*sizeof(double));

	psubSamples = 0;
	pPerturbations = 0;
	PROTECT(rPerturbations = AS_INTEGER(rPerturbations));
	if(!isNull(rPerturbations)) {
		pPerturbations = INTEGER_POINTER(rPerturbations);
		psubSamples = (int*)CATNET_MALLOC(numnodes*numsamples*sizeof(int));
	}

	for(nnode1 = 0; nnode1 < numnodes; nnode1++) {
		psamples1 = pSamples;
		numsamples1 = numsamples;
		if(pPerturbations && psubSamples) {
			numsubsamples = 0;
			for(j = 0; j < numsamples; j++) {
				if(!pPerturbations[j * numnodes + nnode1]) {
					memcpy(psubSamples + numsubsamples*numnodes, pSamples + j*numnodes, numnodes*sizeof(int));
					numsubsamples++;
				}
				psamples1 = psubSamples;
				numsamples1 = numsubsamples;
			}
		}

		for(nnode2 = 0; nnode2 < numnodes; nnode2++) {
	
			memset(pprobs, 0, maxCategories*maxCategories*sizeof(int));
		
			if(nnode2 == nnode1) {
				for(j = 0; j < numsamples1; j++) 
					pprobs[psamples1[j*numnodes + nnode1]]++;
				floglik = 0;
				fsum  = 0;
				for(j = 0; j < pNodeNumCats[nnode1]; j++) {
					fsum += pprobs[j];
					if(pprobs[j] > 0)
						floglik += pprobs[j]*(double)log((double)pprobs[j]);
				}
				if(fsum > 0) {
					floglik -= fsum*(double)log((double)fsum);
					floglik /= fsum;
				}
				klmat[nnode2*numnodes + nnode1] = -floglik;
				continue;
			}

			// estimate logP(nnode1|nnode2)
			for(j = 0; j < numsamples1; j++) 
				pprobs[maxCategories*psamples1[j*numnodes + nnode2] + psamples1[j*numnodes + nnode1]]++;

			floglik = 0;
			fsum  = 0;
			for(i = 0; i < pNodeNumCats[nnode2]; i++) {
				faux = 0;
				for(j = 0; j < pNodeNumCats[nnode1]; j++) {
					faux += pprobs[maxCategories*i+j];
					if(pprobs[maxCategories*i+j] > 0)
						floglik += pprobs[maxCategories*i+j]*(double)log((double)pprobs[maxCategories*i+j]);
				}
				fsum += faux;
				if(faux > 0) {
					floglik -= faux*(double)log((double)faux);
				}
			}
			if(fsum > 0) {
				floglik /= fsum;
			}
			klmat[nnode2*numnodes + nnode1] = -floglik;
		}
	}

	UNPROTECT(2); // rSamples, rPerturbations

	if(psubSamples)
		CATNET_FREE(psubSamples);

	if(pprobs)
		CATNET_FREE(pprobs);

	if(pNodeCats) {
		for(i = 0; i < numnodes; i++) 
			if(pNodeCats[i])
				CATNET_FREE(pNodeCats[i]);
		CATNET_FREE(pNodeCats);
	}

	if(pNodeNumCats) 
		CATNET_FREE(pNodeNumCats);

	if(klmat) {
		PROTECT(rvec = NEW_NUMERIC(numnodes*numnodes));
		pvec = NUMERIC_POINTER(rvec);
		if (pvec && klmat)
			memcpy(pvec, klmat, numnodes*numnodes*sizeof(double));
		UNPROTECT(1);
	}

	CATNET_FREE(klmat);

	//char str[128];
	//sprintf(str, "Mem Balance  %d\n", (int)g_memcounter);
	//printf(str);

	return rvec;

}


SEXP catnetEntropyOrder(SEXP rSamples, SEXP rPerturbations) {

	int *pSamples, *pPerturbations;
	int *psubSamples, numsubsamples, *psamples1, numsamples1;
	int *pNodeNumCats, **pNodeCats, mincat, maxcat, maxCategories, *pprobs;
	int numsamples, numnodes, i, j, k, d, nnode1, nnode2;
	int *porder, *pvec;
	double floglik, faux, fsum, *klmat, *pnoderanks;
	SEXP dim, rvec = R_NilValue;

	if(!isMatrix(rSamples))
		error("Data should be a matrix");
	if(!isNull(rPerturbations) && !isMatrix(rPerturbations))
		error("Perturbations should be a matrix");

	if(isNull(rSamples)) {
		return rvec;
	}
	PROTECT(rSamples = AS_INTEGER(rSamples));
	pSamples = INTEGER(rSamples);
	
	dim = GET_DIM(rSamples);
	numnodes = INTEGER(dim)[0];
	numsamples = INTEGER(dim)[1];

	// pSamples are assumed positive indices
	for(j = 0; j < numnodes*numsamples; j++) {
		pSamples[j]--;
	}

	porder = (int*)CATNET_MALLOC(numnodes*sizeof(int));
	pnoderanks = (double*)CATNET_MALLOC(2*numnodes*sizeof(double));

	// categoies
	pNodeNumCats = (int*)CATNET_MALLOC(numnodes*sizeof(int));
	if (!pNodeNumCats) {
		UNPROTECT(1); //rSamples
		return rvec;
	}
	pNodeCats = (int**)CATNET_MALLOC(numnodes*sizeof(int*));
	if (!pNodeCats) {
		CATNET_FREE(pNodeNumCats);
		UNPROTECT(1); //rSamples
		return rvec;
	}
	memset(pNodeCats, 0, numnodes*sizeof(int*));
	memset(pNodeNumCats, 0, numnodes*sizeof(int));

	maxCategories = 1;
	for(i = 0; i < numnodes; i++) {
		mincat = INT_MAX;
		maxcat = -INT_MAX;
		for(j = 0; j < numsamples; j++) {
			if(pSamples[j*numnodes + i] < mincat)
				mincat = pSamples[j*numnodes + i];
			if(pSamples[j*numnodes + i] > maxcat)
				maxcat = pSamples[j*numnodes + i];
		}
		pNodeNumCats[i] = maxcat - mincat + 1;
		pNodeCats[i] = (int*)CATNET_MALLOC(pNodeNumCats[i]*sizeof(int));
		for(j = 0; j < pNodeNumCats[i]; j++)
			pNodeCats[i][j] = mincat + j;
	}
	for(i = 0; i < numnodes; i++) {
		/* order pNodeNumCats[i] */
		for(j = 0; j < pNodeNumCats[i]; j++) {
			for(k = j + 1; k < pNodeNumCats[i]; k++) {
				if(pNodeCats[i][j] > pNodeCats[i][k]) {
					d = pNodeCats[i][j]; 
					pNodeCats[i][j] = pNodeCats[i][k];
					pNodeCats[i][k] = d;
				}
			}
		} 
		for(j = 0; j < numsamples; j++) {
			for(d = 0; d < pNodeNumCats[i]; d++)
				if(pNodeCats[i][d] == pSamples[j*numnodes + i])
					break;
			pSamples[j*numnodes + i] = d;
		}
		if(maxCategories < pNodeNumCats[i])
			maxCategories = pNodeNumCats[i];
	}

	pprobs = (int*)CATNET_MALLOC(maxCategories*maxCategories*sizeof(int));
	if (!pprobs) {
		CATNET_FREE(pNodeNumCats);
		CATNET_FREE(pNodeCats);
		UNPROTECT(1); //rSamples
		return rvec;
	}

	klmat = (double*)CATNET_MALLOC(numnodes*numnodes*sizeof(double));
	if (!klmat) {
		CATNET_FREE(pNodeNumCats);
		CATNET_FREE(pNodeCats);
		CATNET_FREE(pprobs);
		UNPROTECT(1); //rSamples
		return rvec;
	}

	memset(klmat, 0, numnodes*numnodes*sizeof(double));

	psubSamples = 0;
	pPerturbations = 0;
	PROTECT(rPerturbations = AS_INTEGER(rPerturbations));
	if(!isNull(rPerturbations)) {
		pPerturbations = INTEGER_POINTER(rPerturbations);
		psubSamples = (int*)CATNET_MALLOC(numnodes*numsamples*sizeof(int));
	}

	for(nnode1 = 0; nnode1 < numnodes; nnode1++) {
		psamples1 = pSamples;
		numsamples1 = numsamples;
		if(pPerturbations && psubSamples) {
			numsubsamples = 0;
			for(j = 0; j < numsamples; j++) {
				if(!pPerturbations[j * numnodes + nnode1]) {
					memcpy(psubSamples + numsubsamples*numnodes, pSamples + j*numnodes, numnodes*sizeof(int));
					numsubsamples++;
				}
				psamples1 = psubSamples;
				numsamples1 = numsubsamples;
			}
		}

		for(nnode2 = 0; nnode2 < numnodes; nnode2++) {
	
			memset(pprobs, 0, maxCategories*maxCategories*sizeof(int));
		
			if(nnode2 == nnode1) {
				for(j = 0; j < numsamples1; j++) 
					pprobs[psamples1[j*numnodes + nnode1]]++;
				floglik = 0;
				fsum  = 0;
				for(j = 0; j < pNodeNumCats[nnode1]; j++) {
					fsum += pprobs[j];
					if(pprobs[j] > 0)
						floglik += pprobs[j]*(double)log((double)pprobs[j]);
				}
				if(fsum > 0) {
					floglik -= fsum*(double)log((double)fsum);
					floglik /= fsum;
				}
				klmat[nnode2*numnodes + nnode1] = -floglik;
				continue;
			}

			// estimate logP(nnode1|nnode2)
			for(j = 0; j < numsamples1; j++) 
				pprobs[maxCategories*psamples1[j*numnodes + nnode2] + psamples1[j*numnodes + nnode1]]++;

			floglik = 0;
			fsum  = 0;
			for(i = 0; i < pNodeNumCats[nnode2]; i++) {
				faux = 0;
				for(j = 0; j < pNodeNumCats[nnode1]; j++) {
					faux += pprobs[maxCategories*i+j];
					if(pprobs[maxCategories*i+j] > 0)
						floglik += pprobs[maxCategories*i+j]*(double)log((double)pprobs[maxCategories*i+j]);
				}
				fsum += faux;
				if(faux > 0) {
					floglik -= faux*(double)log((double)faux);
				}
			}
			if(fsum > 0) {
				floglik /= fsum;
			}
			klmat[nnode2*numnodes + nnode1] = -floglik;
		}
	}

	UNPROTECT(2); // rSamples, rPerturbations

	if(psubSamples)
		CATNET_FREE(psubSamples);

	if(pprobs)
		CATNET_FREE(pprobs);

	if(pNodeCats) {
		for(i = 0; i < numnodes; i++) 
			if(pNodeCats[i])
				CATNET_FREE(pNodeCats[i]);
		CATNET_FREE(pNodeCats);
	}

	if(pNodeNumCats) 
		CATNET_FREE(pNodeNumCats);

	for(nnode1 = 0; nnode1 < numnodes; nnode1++) {
		faux = klmat[nnode1*numnodes + nnode1];
		floglik = 0;
		for(nnode2 = 0; nnode2 < numnodes; nnode2++) {
			klmat[nnode2*numnodes + nnode1] = faux - klmat[nnode2*numnodes + nnode1];
			floglik += klmat[nnode2*numnodes + nnode1];
		}
		pnoderanks[nnode1] = floglik;
	}

	_order<double>(pnoderanks, numnodes, porder, 1);

	if (klmat)
		CATNET_FREE(klmat);

	for(i = 0; i < numnodes; i++)
		porder[i]++;

	PROTECT(rvec = NEW_INTEGER(numnodes));
	pvec = INTEGER_POINTER(rvec);
	if (pvec && porder)
		memcpy(pvec, porder, numnodes*sizeof(int));
	UNPROTECT(1);

	CATNET_FREE(porder);
	CATNET_FREE(pnoderanks);

	//char str[128];
	//sprintf(str, "Mem Balance  %d\n", (int)g_memcounter);
	//printf(str);

	return rvec;

}

SEXP catnetKLpairwise(SEXP rSamples, SEXP rPerturbations) {

	int *pSamples, *pPerturbations;
	int *pSamplesPert, numsamplesPert;
	int *pNodeNumCats, **pNodeCats, mincat, maxcat, maxCategories;
	double *pprobs1, *pprobs2;
	int numsamples, numnodes, i, j, k, d, nnode1, nnode2;
	double floglik, faux, fsum, *pvec, *klmat;
	SEXP dim, rvec = R_NilValue;

	if(!isMatrix(rSamples))
		error("Data should be a matrix");
	if(!isNull(rPerturbations) && !isMatrix(rPerturbations))
		error("Perturbations should be a matrix");

	if(isNull(rSamples)) {
		return rvec;
	}
	PROTECT(rSamples = AS_INTEGER(rSamples));
	pSamples = INTEGER(rSamples);
	
	dim = GET_DIM(rSamples);
	numnodes = INTEGER(dim)[0];
	numsamples = INTEGER(dim)[1];	

	if(isNull(rPerturbations)) {
		UNPROTECT(1); //rSamples
		PROTECT(rvec = NEW_NUMERIC(numnodes*numnodes));
		pvec = NUMERIC_POINTER(rvec);
		if (pvec)
			memset(pvec, 0, numnodes*numnodes*sizeof(double));
		UNPROTECT(1);
		return rvec;
	}

	// pSamples are assumed positive indices
	for(j = 0; j < numnodes*numsamples; j++) {
		pSamples[j]--;
	}

	// categoies
	pNodeNumCats = (int*)CATNET_MALLOC(numnodes*sizeof(int));
	if (!pNodeNumCats) {
		UNPROTECT(1); //rSamples
		return rvec;
	}
	pNodeCats = (int**)CATNET_MALLOC(numnodes*sizeof(int*));
	if (!pNodeCats) {
		CATNET_FREE(pNodeNumCats);
		UNPROTECT(1); //rSamples
		return rvec;
	}
	memset(pNodeCats, 0, numnodes*sizeof(int*));
	memset(pNodeNumCats, 0, numnodes*sizeof(int));

	maxCategories = 1;
	for(i = 0; i < numnodes; i++) {
		mincat = INT_MAX;
		maxcat = -INT_MAX;
		for(j = 0; j < numsamples; j++) {
			if(pSamples[j*numnodes + i] < mincat)
				mincat = pSamples[j*numnodes + i];
			if(pSamples[j*numnodes + i] > maxcat)
				maxcat = pSamples[j*numnodes + i];
		}
		pNodeNumCats[i] = maxcat - mincat + 1;
		pNodeCats[i] = (int*)CATNET_MALLOC(pNodeNumCats[i]*sizeof(int));
		for(j = 0; j < pNodeNumCats[i]; j++)
			pNodeCats[i][j] = mincat + j;
	}
	for(i = 0; i < numnodes; i++) {
		/* order pNodeNumCats[i] */
		for(j = 0; j < pNodeNumCats[i]; j++) {
			for(k = j + 1; k < pNodeNumCats[i]; k++) {
				if(pNodeCats[i][j] > pNodeCats[i][k]) {
					d = pNodeCats[i][j]; 
					pNodeCats[i][j] = pNodeCats[i][k];
					pNodeCats[i][k] = d;
				}
			}
		} 
		for(j = 0; j < numsamples; j++) {
			for(d = 0; d < pNodeNumCats[i]; d++)
				if(pNodeCats[i][d] == pSamples[j*numnodes + i])
					break;
			pSamples[j*numnodes + i] = d;
		}
		if(maxCategories < pNodeNumCats[i])
			maxCategories = pNodeNumCats[i];
	}

	pprobs1 = (double*)CATNET_MALLOC(maxCategories*maxCategories*sizeof(double));
	pprobs2 = (double*)CATNET_MALLOC(maxCategories*maxCategories*sizeof(double));

	if (!pprobs1 || !pprobs2) {
		UNPROTECT(1); //rSamples
		CATNET_FREE(pNodeNumCats);
		CATNET_FREE(pNodeCats);
		if (pprobs1)
			CATNET_FREE(pprobs1);
		if (pprobs2)
			CATNET_FREE(pprobs2);
		return rvec;
	}

	pSamplesPert = 0;
	pPerturbations = 0;
	PROTECT(rPerturbations = AS_INTEGER(rPerturbations));
	if(!isNull(rPerturbations)) {
		pPerturbations = INTEGER_POINTER(rPerturbations);
		pSamplesPert = (int*)CATNET_MALLOC(numnodes*numsamples*sizeof(int));
	}

	klmat = (double*)CATNET_MALLOC(numnodes*numnodes*sizeof(double));
	if (!klmat) {
		CATNET_FREE(pNodeNumCats);
		CATNET_FREE(pNodeCats);
		CATNET_FREE(pprobs1);
		CATNET_FREE(pprobs2);
		UNPROTECT(2); // rSamples, rPerturbations
		return rvec;
	}

	memset(klmat, 0, numnodes*numnodes*sizeof(double));

	for(nnode1 = 0; nnode1 < numnodes; nnode1++) {
		numsamplesPert = 0;
		if(pPerturbations && pSamplesPert) {
			for(j = 0; j < numsamples; j++) {
				if(!pPerturbations[j * numnodes + nnode1]) {
					memcpy(pSamplesPert + numsamplesPert*numnodes, pSamples + j*numnodes, numnodes*sizeof(int));
					numsamplesPert++;
				}
			}
		}
		//printf("\nnnode = %d (%d)\n", nnode1, numsamplesPert);
		for(nnode2 = 0; nnode2 < numnodes; nnode2++) {
	
			if(nnode1 == nnode2)
				continue;

			memset(pprobs2, 0, maxCategories*maxCategories*sizeof(double));
			// estimate logP(nnode1|nnode2) for the whole sample
			for(j = 0; j < numsamples; j++) 
				pprobs2[maxCategories*pSamples[j*numnodes + nnode2] + pSamples[j*numnodes + nnode1]]+=1;

			memset(pprobs1, 0, maxCategories*maxCategories*sizeof(double));
			// estimate logP(nnode1|nnode2) for the perturbed sub-sample only
			for(j = 0; j < numsamplesPert; j++) 
				pprobs1[maxCategories*pSamplesPert[j*numnodes + nnode2] + pSamplesPert[j*numnodes + nnode1]]+=1;

			for(i = 0; i < pNodeNumCats[nnode2]; i++) {
				fsum = 0;
				for(j = 0; j < pNodeNumCats[nnode1]; j++)
					fsum += pprobs1[maxCategories*i+j];
				if(fsum <= 0)
					continue;
				faux = 1 / fsum;
				for(j = 0; j < pNodeNumCats[nnode1]; j++)
					pprobs1[maxCategories*i+j] *= faux;
			}
			for(i = 0; i < pNodeNumCats[nnode2]; i++) {
				fsum = 0;
				for(j = 0; j < pNodeNumCats[nnode1]; j++)
					fsum += pprobs2[maxCategories*i+j];
				if(fsum <= 0)
					continue;
				faux = 1 / fsum;
				for(j = 0; j < pNodeNumCats[nnode1]; j++)
					pprobs2[maxCategories*i+j] *= faux;
			}

			floglik = 0;
			for(i = 0; i < pNodeNumCats[nnode2]; i++) {
				for(j = 0; j < pNodeNumCats[nnode1]; j++) {
					if(pprobs1[maxCategories*i+j] > 0 && pprobs2[maxCategories*i+j] > 0)
						floglik += pprobs1[maxCategories*i+j]*
						(double)log((double)pprobs1[maxCategories*i+j] / (double)pprobs2[maxCategories*i+j]);
					else if(pprobs1[maxCategories*i+j] != 0 && pprobs2[maxCategories*i+j] == 0)
						floglik = FLT_MAX;
				}
			}
			klmat[nnode2*numnodes + nnode1] += floglik;
		}
	}

	UNPROTECT(2); // rSamples, rPerturbations

	if(pSamplesPert)
		CATNET_FREE(pSamplesPert);

	if(pprobs1)
		CATNET_FREE(pprobs1);
	if(pprobs2)
		CATNET_FREE(pprobs2);

	if(pNodeCats) {
		for(i = 0; i < numnodes; i++) 
			if(pNodeCats[i])
				CATNET_FREE(pNodeCats[i]);
		CATNET_FREE(pNodeCats);
	}

	if(pNodeNumCats) 
		CATNET_FREE(pNodeNumCats);

	PROTECT(rvec = NEW_NUMERIC(numnodes*numnodes));
	pvec = NUMERIC_POINTER(rvec);
	if (pvec && klmat)
		memcpy(pvec, klmat, numnodes*numnodes*sizeof(double));
	UNPROTECT(1);

	if (klmat)
		CATNET_FREE(klmat);

	//char str[128];
	//sprintf(str, "Mem Balance  %d\n", (int)g_memcounter);
	//printf(str);

	return rvec;

}

SEXP catnetPearsonPairwise(SEXP rSamples, SEXP rPerturbations) {

	int *pSamples, *pPerturbations;
	int *pSamplesPert, numsamplesPert;
	int *pNodeNumCats, **pNodeCats, mincat, maxcat, maxCategories;
	double *pprobs1, *pprobs2;
	int numsamples, numnodes, i, j, k, d, nnode1, nnode2;
	double floglik, faux, fsum, *pvec, *klmat;
	SEXP dim, rvec = R_NilValue;

	if(!isMatrix(rSamples))
		error("Data should be a matrix");
	if(!isNull(rPerturbations) && !isMatrix(rPerturbations))
		error("Perturbations should be a matrix");

	if(isNull(rSamples)) {
		return rvec;
	}
	PROTECT(rSamples = AS_INTEGER(rSamples));
	pSamples = INTEGER(rSamples);
	
	dim = GET_DIM(rSamples);
	numnodes = INTEGER(dim)[0];
	numsamples = INTEGER(dim)[1];	

	if(isNull(rPerturbations)) {
		UNPROTECT(1); //rSamples
		PROTECT(rvec = NEW_NUMERIC(numnodes*numnodes));
		pvec = NUMERIC_POINTER(rvec);
		if (pvec)
			memset(pvec, 0, numnodes*numnodes*sizeof(double));
		UNPROTECT(1);
		return rvec;
	}

	// pSamples are assumed positive indices
	for(j = 0; j < numnodes*numsamples; j++) {
		pSamples[j]--;
	}

	// categoies
	pNodeNumCats = (int*)CATNET_MALLOC(numnodes*sizeof(int));
	if (!pNodeNumCats) {
		UNPROTECT(1); //rSamples
		return rvec;
	}
	pNodeCats = (int**)CATNET_MALLOC(numnodes*sizeof(int*));
	if (!pNodeCats) {
		CATNET_FREE(pNodeNumCats);
		UNPROTECT(1); //rSamples
		return rvec;
	}
	memset(pNodeCats, 0, numnodes*sizeof(int*));
	memset(pNodeNumCats, 0, numnodes*sizeof(int));

	maxCategories = 1;
	for(i = 0; i < numnodes; i++) {
		mincat = INT_MAX;
		maxcat = -INT_MAX;
		for(j = 0; j < numsamples; j++) {
			if(pSamples[j*numnodes + i] < mincat)
				mincat = pSamples[j*numnodes + i];
			if(pSamples[j*numnodes + i] > maxcat)
				maxcat = pSamples[j*numnodes + i];
		}
		pNodeNumCats[i] = maxcat - mincat + 1;
		pNodeCats[i] = (int*)CATNET_MALLOC(pNodeNumCats[i]*sizeof(int));
		for(j = 0; j < pNodeNumCats[i]; j++)
			pNodeCats[i][j] = mincat + j;
	}
	for(i = 0; i < numnodes; i++) {
		/* order pNodeNumCats[i] */
		for(j = 0; j < pNodeNumCats[i]; j++) {
			for(k = j + 1; k < pNodeNumCats[i]; k++) {
				if(pNodeCats[i][j] > pNodeCats[i][k]) {
					d = pNodeCats[i][j]; 
					pNodeCats[i][j] = pNodeCats[i][k];
					pNodeCats[i][k] = d;
				}
			}
		} 
		for(j = 0; j < numsamples; j++) {
			for(d = 0; d < pNodeNumCats[i]; d++)
				if(pNodeCats[i][d] == pSamples[j*numnodes + i])
					break;
			pSamples[j*numnodes + i] = d;
		}
		if(maxCategories < pNodeNumCats[i])
			maxCategories = pNodeNumCats[i];
	}

	pprobs1 = (double*)CATNET_MALLOC(maxCategories*maxCategories*sizeof(double));
	pprobs2 = (double*)CATNET_MALLOC(maxCategories*maxCategories*sizeof(double));
	if (!pprobs1 || !pprobs2) {
		UNPROTECT(1); //rSamples
		CATNET_FREE(pNodeNumCats);
		CATNET_FREE(pNodeCats);
		if (pprobs1)
			CATNET_FREE(pprobs1);
		if (pprobs2)
			CATNET_FREE(pprobs2);
		return rvec;
	}

	pSamplesPert = 0;
	pPerturbations = 0;
	PROTECT(rPerturbations = AS_INTEGER(rPerturbations));
	if(!isNull(rPerturbations)) {
		pPerturbations = INTEGER_POINTER(rPerturbations);
		pSamplesPert = (int*)CATNET_MALLOC(numnodes*numsamples*sizeof(int));
	}

	klmat = (double*)CATNET_MALLOC(numnodes*numnodes*sizeof(double));
	if (!klmat) {
		CATNET_FREE(pNodeNumCats);
		CATNET_FREE(pNodeCats);
		CATNET_FREE(pprobs1);
		CATNET_FREE(pprobs2);
		UNPROTECT(2); //rSamples, rPerturbations
		return rvec;
	}
	
	memset(klmat, 0, numnodes*numnodes*sizeof(double));

	for(nnode1 = 0; nnode1 < numnodes; nnode1++) {
		numsamplesPert = 0;
		if(pPerturbations && pSamplesPert) {
			for(j = 0; j < numsamples; j++) {
				if(!pPerturbations[j * numnodes + nnode1]) {
					memcpy(pSamplesPert + numsamplesPert*numnodes, pSamples + j*numnodes, numnodes*sizeof(int));
					numsamplesPert++;
				}
			}
		}
		//printf("\nnnode = %d (%d)\n", nnode1, numsamplesPert);
		for(nnode2 = 0; nnode2 < numnodes; nnode2++) {
	
			if(nnode1 == nnode2)
				continue;

			memset(pprobs2, 0, maxCategories*maxCategories*sizeof(double));
			// estimate logP(nnode1|nnode2) for the whole sample
			for(j = 0; j < numsamples; j++) 
				pprobs2[maxCategories*pSamples[j*numnodes + nnode2] + pSamples[j*numnodes + nnode1]]+=1;

			memset(pprobs1, 0, maxCategories*maxCategories*sizeof(double));
			// estimate logP(nnode1|nnode2) for the perturbed sub-sample only
			for(j = 0; j < numsamplesPert; j++) 
				pprobs1[maxCategories*pSamplesPert[j*numnodes + nnode2] + pSamplesPert[j*numnodes + nnode1]]+=1;

			for(i = 0; i < pNodeNumCats[nnode2]; i++) {
				fsum = 0;
				for(j = 0; j < pNodeNumCats[nnode1]; j++)
					fsum += pprobs2[maxCategories*i+j];
				if(fsum <= 0)
					continue;
				faux = 1 / fsum;
				for(j = 0; j < pNodeNumCats[nnode1]; j++)
					pprobs2[maxCategories*i+j] *= faux;
			}

			floglik = 0;
			for(i = 0; i < pNodeNumCats[nnode2]; i++) {
				fsum = 0;
				for(j = 0; j < pNodeNumCats[nnode1]; j++)
					fsum += pprobs1[maxCategories*i+j];
				if(fsum <= 0)
					continue;
				for(j = 0; j < pNodeNumCats[nnode1]; j++) {
					faux = pprobs1[maxCategories*i+j] - fsum*pprobs2[maxCategories*i+j];
					if(pprobs2[maxCategories*i+j] > 0)
						floglik += (double)(faux*faux) / (double)(fsum*pprobs2[maxCategories*i+j]);
					else if(faux != 0 && pprobs2[maxCategories*i+j] == 0)
						floglik = FLT_MAX;
				}
			}
			klmat[nnode2*numnodes + nnode1] += floglik;
		}
	}

	UNPROTECT(2); //rSamples, rPerturbations

	if(pSamplesPert)
		CATNET_FREE(pSamplesPert);

	if(pprobs1)
		CATNET_FREE(pprobs1);
	if(pprobs2)
		CATNET_FREE(pprobs2);

	if(pNodeCats) {
		for(i = 0; i < numnodes; i++) 
			if(pNodeCats[i])
				CATNET_FREE(pNodeCats[i]);
		CATNET_FREE(pNodeCats);
	}

	if(pNodeNumCats) 
		CATNET_FREE(pNodeNumCats);

	PROTECT(rvec = NEW_NUMERIC(numnodes*numnodes));
	pvec = NUMERIC_POINTER(rvec);
	if (pvec && klmat)
		memcpy(pvec, klmat, numnodes*numnodes*sizeof(double));
	UNPROTECT(1);

	if (klmat)
		CATNET_FREE(klmat);

	//char str[128];
	//sprintf(str, "Mem Balance  %d\n", (int)g_memcounter);
	//printf(str);

	return rvec;

}

} // extern "C"


double *catnetPairwiseCondLikelihood(SEXP rSamples, SEXP rPerturbations) {

	int *pSamples, *pRSamples, *pPerturbations;
	int *pSamplesPert, numsamplesPert;
	int *pNodeNumCats, **pNodeCats, mincat, maxcat, maxCategories;
	double *pprobs;
	int numsamples, numnodes, i, j, k, d, nnode1, nnode2, ncount;
	double floglik, fsum, *matPairs, fmin, fmax;
	SEXP dim;

	if(!isMatrix(rSamples))
		error("Data should be a matrix");
	if(!isNull(rPerturbations) && !isMatrix(rPerturbations))
		error("Perturbations should be a matrix");

	if(isNull(rSamples)) {
		return 0;
	}
	PROTECT(rSamples = AS_INTEGER(rSamples));
	pRSamples = INTEGER(rSamples);
	
	dim = GET_DIM(rSamples);
	numnodes = INTEGER(dim)[0];
	numsamples = INTEGER(dim)[1];	

	// pSamples are assumed positive indices
	pSamples = (int*)CATNET_MALLOC(numnodes*numsamples*sizeof(int));
	if (!pSamples) {
		UNPROTECT(1); //rSamples
		return 0;
	}
	memcpy(pSamples, pRSamples, numnodes*numsamples*sizeof(int));
	UNPROTECT(1); //rSamples

	for(j = 0; j < numnodes*numsamples; j++) {
		pSamples[j]--;
	}

	// categoies
	pNodeNumCats = (int*)CATNET_MALLOC(numnodes*sizeof(int));
	if (!pNodeNumCats) {
		CATNET_FREE(pSamples);
		return 0;
	}
	pNodeCats = (int**)CATNET_MALLOC(numnodes*sizeof(int*));
	if (!pNodeCats) {
		CATNET_FREE(pSamples);
		CATNET_FREE(pNodeNumCats);
		return 0;
	}
	memset(pNodeCats, 0, numnodes*sizeof(int*));
	memset(pNodeNumCats, 0, numnodes*sizeof(int));

	maxCategories = 1;
	for(i = 0; i < numnodes; i++) {
		mincat = INT_MAX;
		maxcat = -INT_MAX;
		for(j = 0; j < numsamples; j++) {
			if(pSamples[j*numnodes + i] < mincat)
				mincat = pSamples[j*numnodes + i];
			if(pSamples[j*numnodes + i] > maxcat)
				maxcat = pSamples[j*numnodes + i];
		}
		pNodeNumCats[i] = maxcat - mincat + 1;
		pNodeCats[i] = (int*)CATNET_MALLOC(pNodeNumCats[i]*sizeof(int));
		for(j = 0; j < pNodeNumCats[i]; j++)
			pNodeCats[i][j] = mincat + j;
	}
	for(i = 0; i < numnodes; i++) {
		/* order pNodeNumCats[i] */
		for(j = 0; j < pNodeNumCats[i]; j++) {
			for(k = j + 1; k < pNodeNumCats[i]; k++) {
				if(pNodeCats[i][j] > pNodeCats[i][k]) {
					d = pNodeCats[i][j]; 
					pNodeCats[i][j] = pNodeCats[i][k];
					pNodeCats[i][k] = d;
				}
			}
		} 
		for(j = 0; j < numsamples; j++) {
			for(d = 0; d < pNodeNumCats[i]; d++)
				if(pNodeCats[i][d] == pSamples[j*numnodes + i])
					break;
			pSamples[j*numnodes + i] = d;
		}
		if(maxCategories < pNodeNumCats[i])
			maxCategories = pNodeNumCats[i];
	}

	pprobs = (double*)CATNET_MALLOC(maxCategories*maxCategories*sizeof(double));
	if (!pprobs) {
		CATNET_FREE(pSamples);
		CATNET_FREE(pNodeNumCats);
		CATNET_FREE(pNodeCats);
		return 0;
	}

	matPairs = (double*)CATNET_MALLOC(numnodes*numnodes*sizeof(double));
	if (!matPairs) {
		CATNET_FREE(pSamples);
		CATNET_FREE(pNodeNumCats);
		CATNET_FREE(pNodeCats);
		CATNET_FREE(pprobs);
		return 0;	
	}

	memset(matPairs, 0, numnodes*numnodes*sizeof(double));

	pSamplesPert = 0;
	pPerturbations = 0;
	PROTECT(rPerturbations = AS_INTEGER(rPerturbations));
	if(!isNull(rPerturbations)) {
		pPerturbations = INTEGER_POINTER(rPerturbations);
		pSamplesPert = (int*)CATNET_MALLOC(numnodes*numsamples*sizeof(int));
	}

	for(nnode1 = 0; nnode1 < numnodes; nnode1++) {
		numsamplesPert = 0;
		if(pPerturbations && pSamplesPert) {
			for(j = 0; j < numsamples; j++) {
				if(!pPerturbations[j * numnodes + nnode1]) {
					memcpy(pSamplesPert + numsamplesPert*numnodes, pSamples + j*numnodes, numnodes*sizeof(int));
					numsamplesPert++;
				}
			}
		}

		for(nnode2 = 0; nnode2 < numnodes; nnode2++) {
	
			if(nnode1 == nnode2)
				continue;

			ncount = 0;
			memset(pprobs, 0, maxCategories*maxCategories*sizeof(double));
			// estimate logP(nnode1|nnode2)
			if(pPerturbations && pSamplesPert) {
				for(j = 0; j < numsamplesPert; j++) {
					pprobs[maxCategories*pSamplesPert[j*numnodes + nnode2] + pSamplesPert[j*numnodes + nnode1]]++; 
					ncount++;
				}
			} 
			else {
				for(j = 0; j < numsamples; j++) {
					pprobs[maxCategories*pSamples[j*numnodes + nnode2] + pSamples[j*numnodes + nnode1]]++; 
					ncount++;
				}
			}

			floglik = 0;
			for(i = 0; i < pNodeNumCats[nnode2]; i++) {
				fsum = 0;
				fmin = 0;
				for(j = 0; j < pNodeNumCats[nnode1]; j++) {
					fsum += pprobs[maxCategories*i+j];
					if(pprobs[maxCategories*i+j] > 0)
						fmin += (double)pprobs[maxCategories*i+j] * (double)log((double)pprobs[maxCategories*i+j]);
				}
				floglik += fmin;
				if(fsum > 0)
					floglik -= fsum*log(fsum);
			}
			if(ncount > 1 && floglik > (double)-FLT_MAX)
				floglik /= (double)ncount;
			matPairs[nnode1*numnodes + nnode2] = floglik;
		}
		fsum = 0; fmin = FLT_MAX; fmax = -FLT_MAX;
		for(nnode2 = 0; nnode2 < numnodes; nnode2++) {
			fsum += matPairs[nnode1*numnodes + nnode2];
			if(fmin > matPairs[nnode1*numnodes + nnode2])
				fmin = matPairs[nnode1*numnodes + nnode2];
			if(fmax < matPairs[nnode1*numnodes + nnode2])
				fmax = matPairs[nnode1*numnodes + nnode2];
		}

		fsum = 1;
		if(fmax-fmin>0)
			fsum = 1 / (fmax-fmin);
		for(nnode2 = 0; nnode2 < numnodes; nnode2++) 
			 matPairs[nnode1*numnodes + nnode2] = (matPairs[nnode1*numnodes + nnode2] - fmin)*fsum; 
	}

	UNPROTECT(1); //rPerturbations

	if(pSamplesPert)
		CATNET_FREE(pSamplesPert);

	if(pSamples)
		CATNET_FREE(pSamples);

	if(pprobs)
		CATNET_FREE(pprobs);

	if(pNodeCats) {
		for(i = 0; i < numnodes; i++) 
			if(pNodeCats[i])
				CATNET_FREE(pNodeCats[i]);
		CATNET_FREE(pNodeCats);
	}

	if(pNodeNumCats) 
		CATNET_FREE(pNodeNumCats);

	//char str[128];
	//sprintf(str, "Mem Balance  %d\n", (int)g_memcounter);
	//printf(str);

	return matPairs;

}
