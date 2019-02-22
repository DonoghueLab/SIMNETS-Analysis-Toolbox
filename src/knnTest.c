#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "knnHelpers.h"
#include "randperm.h"

size_t n = 9;
size_t n_cols = 3;
size_t k = 3;
size_t n_perm = 100000;
/*
double data[] = {	1.0,	1.2,	0.9,	3.4,	3.6,	2.9,	2.2,	2.1,	2.5,
					2.0,	2.2,	2.1,	0.1,	-0.2,	0.3,	1.5,	1.8,	1.2,
					0.6,	1.4,	0.9,	2.1,	0.2,	1.6,	1.2,	0.8,	2.3};
//*/
double data[] = {	1.0,	2.0,	0.6,
					1.2,	2.2,	1.4,
					0.9,	2.1,	0.9,
					3.4,	0.1,	2.1,
					3.6,	-0.2,	0.2,
					2.9,	0.3,	1.6,
					2.2,	1.5,	1.2,
					2.1,	1.8,	0.8,
					2.5,	1.2,	2.3};
double dmat[] = {	1.1,	0.9,	2.1,	2.5,	4.8,	3.4,	2.1,	3.7,
												0.5,	1.9,	2.9,	4.2,	3.1,	4.0,	3.7,
															2.7,	1.8,	3.3,	5.2,	2.8,	3.7,
																		0.3,	2.2,	1.9,	0.5,	3.7,
																					3.7,	2.3,	0.2,	3.7,
																								0.4,	3.8,	1.3,
																											2.9,	0.8,
																														3.7};
long long classes[] =  {1, 1, 1, 3, 3, 3, 7, 7, 7};

int main (int argc, char const *argv[])
{
//	time_t t;
	double perc, *pperm, *vm;
	long long *knnclass, *uclasses;
	pperm = malloc(n_perm * sizeof(double));
	knnclass = malloc(n * sizeof(long long));
	uclasses = malloc(n * sizeof(long long));

	ssize_t n_unique_classes = findUniqueItemsLL(classes, n, uclasses);
	vm = calloc(n * n_unique_classes, sizeof(double));

	double *dist_mat_tril = malloc(n * (n - 1) / 2 * sizeof(double));
	Rpdist(data, n, n_cols, dist_mat_tril);


	srand((unsigned) time(NULL));
	perc = doKnnClassification(dist_mat_tril, n, classes, k, n_perm, knnclass, vm, pperm);
	printf("\nClassification was correct %4.1f%%.\n", perc);

	if (n_perm > 0) {
		size_t i_up = floor(n_perm * 0.975),
			i_low = floor(n_perm * 0.025),
			i_med = floor(n_perm * 0.5);
		printf("%zu %zu %zu\n", i_low, i_med, i_up);
		printf("Permutation correct: [%5.1f - %5.1f - %5.1f]\n", pperm[i_low], pperm[i_med], pperm[i_up]);
	}


	if (knnclass) {
		for (size_t i = 0; i < n; ++i) {
			printf("\n%3zu: %3lld -> %3lld\n", i, classes[i], knnclass[i]);
			for(size_t j = 0; j < n_unique_classes; ++j)
			{
				printf("%5.2f ", vm[i*n_unique_classes + j]);
			}
		}
	}

	size_t m = 6;
	printf("\n\n");
	for(size_t row = 0; row < m; ++row)
	{
		printf("Row: %zu\t", row);
		for(size_t col = 0; col < m; ++col)
		{
			if (col == row) {
				printf("   -");
			} else {
				printf("%4zu", sub2ind_tril(m, row, col));
			}
		}
		printf("\n");
	}
	/*
	size_t n = 1001, perms = 100000;
	size_t * perm_ind = malloc(n * sizeof(size_t));
	double first = 0, last = 0;
	for(size_t j = 0; j < perms; ++j) {
		randperm(n, perm_ind);
/*
		printf("1: ");
		for (size_t i = 0; i < n; ++i) {
			printf("%3zu ", perm_ind[i]);
		}
		printf("\n");
//* /
		first += perm_ind[0];
		last += perm_ind[n-1];
	}

	first = first / perms;
	last = last / perms;
	printf("Avg first: %f, last: %f.\n", first, last);
	free(perm_ind);

//*/

	free(uclasses);
	free(knnclass);
	free(pperm);
	return 0;
}
