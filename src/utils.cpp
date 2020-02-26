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
 * utils.c
 *
 *  Created on: Nov 16, 2009
 *      Author: Nikolay Balov
 */

/* 
 * version 1.15.1  12dec2016
 */

#include "utils.h"

using namespace std;

extern size_t g_memcounter;

void * CATNET_MALLOC(size_t nsize) { 
	if(nsize <= 0)
		return 0;
	void *pMem;
	pMem = malloc(nsize);
	if(!pMem) {
		error("Insufficient memory");
		return 0;
	}
	return pMem;
	// memmonitor
	g_memcounter += nsize;
	pMem = malloc(sizeof(int) + nsize);
	if(!pMem) {
		error("Insufficient memory");
		return 0;
	}
	*(int*)pMem = nsize;
	pMem = (void*)((int*)pMem + 1);
	//char str[128];
	//sprintf(str, "+%d    %d        %p\n", (int)nsize, (int)g_memcounter, pMem);
	//fprintf(g_hf,str);
	//printf(str);
	return pMem;
}

void CATNET_FREE(void *pMem) {
	if(!pMem)	
		return;
	free(pMem);
	return;
	// memmonitor
	pMem = (void*)((int*)pMem-1);
	size_t nsize = *((int*)pMem);
	g_memcounter -= nsize;
	//char str[128];
	//sprintf(str, "-%d    %d\n", (int)nsize, (int)g_memcounter);
	//printf(str);
	free(pMem);
}

void CATNET_MEM_ERR() {
	// generate R-errorx
	error("Insufficient memory");
}
