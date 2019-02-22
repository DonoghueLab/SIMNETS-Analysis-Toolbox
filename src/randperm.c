/*==========================================================
    random array helpers

	copyright notice for randperm():
		Copyright (c) 1990 Michael E. Hohmeyer,
				hohmeyer@icemcfd.com
		Permission is granted to modify and re-distribute this code in any manner
		as long as this notice is preserved.  All standard disclaimers apply.

    Jonas Zimmermann (c) 2013

    @author Jonas Zimmermann
    Copyright (c) Jonas Zimmermann, Brown University. All rights reserved.

 *========================================================*/

#include "randperm.h"
#include <stdlib.h>
#include <time.h>

/* Swap elements I and J in array V.  */
static void
swap (size_t *v, size_t i, size_t j)
{
  size_t t = v[i];
  v[i] = v[j];
  v[j] = t;
}


void randperm(size_t n, size_t perm[])
{
	for (size_t i = 0; i < n; i++) {
		perm[i] = i;
	}
	for (size_t i = 0; i < n; i++) {
		size_t j = rand() % (n - i) + i;
		swap(perm, i, j);
	}
}

void randsamp(size_t n, size_t k, size_t samp[])
{
	for (size_t i = 0; i < k; ++i) {
		samp[i] = rand()%n;
	}
}
