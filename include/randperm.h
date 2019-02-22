/*==========================================================
    random array helpers

    Jonas Zimmermann (c) 2013

    @author Jonas Zimmermann
    Copyright (c) Jonas Zimmermann, Brown University. All rights reserved.

 *========================================================*/
#ifndef RANDPERM_H_K0NTLYJZ
#define RANDPERM_H_K0NTLYJZ


#include <stdlib.h>

void randperm(size_t n, size_t perm[]);
/*
Calculates a random permutaion of n integers [0 .. n-1]
Input args:		n			number of items
Output args:	perm		preallocated n-array containing the permutation
*/

void randsamp(size_t n, size_t k, size_t samp[]);
/*
Calculates a random sample drawn k times from n numbers [0 .. n-1], with replacement.
*/
#endif /* end of include guard: RANDPERM_H_K0NTLYJZ */
