#ifdef MATLAB_MEX_FILE
#include "mex.h"
#define malloc mxMalloc
#define calloc mxCalloc
#define realloc mxRealloc
#define free mxFree
#define printf mexPrintf
#else
#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <stdlib.h>

#define mexMakeMemoryPersistent(a)
#define mexEvalString(a)
#endif

#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/types.h>
#include "knnHelpers.h"
#include "randperm.h"

typedef struct
{
	double d;
	size_t i;
} double_with_index_t;

#define FREE_IF_NOT_NULL(p) 	if ( p ) {\
		free ( p ); \
		p = NULL; \
	};


static size_t knnFind(const long long *classes, const size_t n_classes, const long long *nn_classes, const double *weights, const size_t n_nn, double *votes);
static int cmpfunc (const void * a, const void * b);
static int dcmpfunc (const void * a, const void * b);
static int compare_doubles_with_index(const void *a, const void *b);
static double dist_c(const double* x, const double *y, const size_t p);


size_t sub2ind_tril(size_t n_rows, size_t row, size_t col)
{
	if (col > row) {
		size_t tmp = row;
		row = col;
		col = tmp;
	}
	return row - 1 + (col ) * (n_rows-1) - col * (col + 1) / 2;
}

static size_t knnFind(const long long *classes, const size_t n_classes, const long long *nn_classes, const double *weights, const size_t n_nn, double *votes)
{
	size_t max_i = 0;
	double max_v = 0;

	for (size_t i = 0; i < n_classes; ++i) {
		votes[i] = 0;
		for (size_t j = 0; j < n_nn; ++j) {
			if (nn_classes[j] == classes[i]) {
				votes[i] += weights[j];
			}
		}
		if (votes[i] > max_v) {
			max_v = votes[i];
			max_i = i;
		}
	}

	return max_i;
}

static int cmpfunc (const void * a, const void * b)
{
	const long long *p_a = a;
	const long long *p_b = b;

	if (*p_a < *p_b) { return -1;}
	if (*p_a > *p_b) { return 1;}
	return 0;
}

static int dcmpfunc (const void * a, const void * b)
{
	const double *p_a = a;
	const double *p_b = b;

	return (*p_a > *p_b) - (*p_a < *p_b);
}

static int compare_doubles_with_index(const void *a, const void *b)
{
	const double_with_index_t *da = (const double_with_index_t *) a;
	const double_with_index_t *db = (const double_with_index_t *) b;

	return (da->d > db->d) - (da->d < db->d);
}

static double dist(const double* x, const double *y, const size_t nx, const size_t ny, const size_t p)
{
	double distance = 0, xi, yi, diff;
#ifdef PRINT_DEBUG
	printf("Item 1:\t");
	for (size_t i = 0; i < p; i++) {
		printf("%4.1f\t", x[i*nx]);
	}
	printf("\nItem 2:\t");
	for (size_t i = 0; i < p; i++) {
		printf("%4.1f\t", y[i*ny]);
	}
	printf("\n\n");
#endif
	for (size_t i = 0; i < p; i++) {
		//i*nx is the ith feature of the nxth row
		xi = x[i*nx];
		yi = y[i*ny];
		if (isnan(xi) || isnan(yi)) {
			continue;
		}
		diff = xi - yi;
		distance += (diff*diff);
	}
	return sqrt(distance);
}

static double dist_c(const double* x, const double *y, const size_t p)
{
	double distance = 0, xi, yi, diff;
	for (size_t i = 0; i < p; i++) {
		xi = x[i];
		yi = y[i];
		if (isnan(xi) || isnan(yi)) {
			continue;
		}
		diff = xi - yi;
		distance += (diff*diff);
	}
	return sqrt(distance);
}


// C implementation of pdist
// X is a matrix represented as a vector.  Columns of X are features
// n_rows is the number of observations (rows) in X
// n_cols is the length of the vectors in X
void Rpdist(const double* X, const size_t n_rows, const size_t n_cols, double* distances)
{
	for (size_t i = 0; i < n_rows - 1; i++) {
		const double *x = X + i;  //a pointer to the beginning of the ith row of X
		for (size_t j = i + 1; j < n_rows; j++) {
			const double *y = X + j;
			distances[sub2ind_tril(n_rows, i, j)] = dist(x, y, n_rows, n_rows, n_cols);
		}
	}
}

// C implementation of pdist
// X is a matrix represented as a vector.  Rows of X are features
// n_rows is the length of the vectors in X
// n_cols is the number of observations (rows) in X
void Rpdist_c(const double* X, const size_t n_rows, const size_t n_cols, double* distances)
{
	for (size_t i = 0; i < n_cols - 1; i++) {
		const double *x = X + i * n_rows;  //a pointer to the beginning of the ith row of X
		for (size_t j = i + 1; j < n_cols; j++) {
			const double *y = X + j * n_rows;
			distances[sub2ind_tril(n_cols, i, j)] = dist_c(x, y, n_rows);
		}
	}
}


ssize_t findUniqueItemsLL(const long long *list, const size_t n, long long *ulist)
{
	if (ulist) {
		memcpy(ulist, list, n * sizeof(long long));
		qsort(ulist, n, sizeof(long long), cmpfunc);
		ssize_t n_unique_classes = 1;
		for (size_t i=1; i < n; i++) {
			if (ulist[i-1]^ulist[i]) {
				ulist[n_unique_classes ++] = ulist[i];
			}
		}
		return n_unique_classes;
	}
	return -1;
}


double doKnnClassification(const double *dist_mat_tril, const size_t n, const long long *true_class, const size_t k, const size_t n_perm, long long *knn_class, double *vote_matrix, double *per_perm)
{
	long long *uclasses = malloc(n * sizeof(long long));
	size_t *_knn_class = malloc(n * sizeof(size_t));
	long long *nnclasses = malloc((n - 1) * sizeof(long long));

	ssize_t n_unique_classes = findUniqueItemsLL(true_class, n, uclasses);

	double *v_m = calloc(n_unique_classes, sizeof(double));

	double *p_perm = NULL;
	size_t *n_corr_perm = NULL;
	size_t *perm_ind = NULL;
	if (n_perm > 0) {
		p_perm = calloc(n_perm, sizeof(double));
		n_corr_perm = calloc(n_perm, sizeof(size_t));
		perm_ind = calloc(n, sizeof(size_t));
	}
	double_with_index_t *d_line = malloc((n - 1) * sizeof(double_with_index_t));
	double *weights = malloc((n - 1) * sizeof(double));
	double num_correct = 0;
	size_t ind, this_k;

	time_t reference_time;
	time(&reference_time);

	for (size_t i = 0; i < n; ++i) {
		// load current item's distances
		for (size_t j = 0; j < i; ++j) {
			d_line[j].d = dist_mat_tril[sub2ind_tril(n, j, i)];
			d_line[j].i = j;
		}
		for (size_t j = 0; j < n - i - 1 ; ++j) {
			d_line[j + i].d = dist_mat_tril[sub2ind_tril(n, j + i + 1, i)];
			d_line[j + i].i = i + j + 1;
		}

		qsort(d_line, n - 1, sizeof(*d_line), compare_doubles_with_index);

		this_k = k;
		for (size_t j = 0; j < this_k; ++j) {
			ind = d_line[j].i;
			nnclasses[j] = true_class[ind];
			weights[j] = 1.0 / d_line[j].d;
			// If we're at the k-th item and there are still more items, check if they are at the same distance. If so,
			// also add the to neighbour's list
			if ((j == this_k - 1) && (j < n - 2) && (d_line[j].d == d_line[j + 1].d)) {
				++ this_k;
			}
		}

		_knn_class[i] = knnFind(uclasses, n_unique_classes, nnclasses, weights, this_k, v_m);
		if (uclasses[_knn_class[i]] == true_class[i]) {
			++ num_correct;
		}
		if (vote_matrix) {
			memcpy(vote_matrix + i*n_unique_classes, v_m, n_unique_classes * sizeof(double));
		}

		// now do permutation
		if (n_perm > 0) {
			size_t this_class;
			// for each object to be tested we initialize the random number generator with the same time determined
			// at start of this function. This ensures that every item sees the same permutations of true classes,
			// which should ensure the same outcome as if repeatedly calling KNN calssify and permuting true_class.
			srand((unsigned) reference_time);
			for (size_t j = 0; j < n_perm; ++j) {
				randperm(n, perm_ind);

				for (size_t kk = 0; kk < this_k; ++kk) {
					nnclasses[kk] = true_class[perm_ind[d_line[kk].i]];
				}
				this_class = knnFind(uclasses, n_unique_classes, nnclasses, weights, this_k, v_m);
				if (uclasses[this_class] == true_class[perm_ind[i]]) {
					++ n_corr_perm[j];
				}
			}
		}

#ifdef PRINT_DEBUG
		printf("\nItem i: %4zu ", i);
		for (size_t j = 0; j < n-1; ++j) {
			printf("%5.1f (%2zu |%2zi) ", d_line[j].d, d_line[j].i, true_class[d_line[j].i]);
		}
		printf("\tClass: %zi (%zu)\tVotes: ", uclasses[_knn_class[i]], _knn_class[i]);
		for (size_t j = 0; j < n_unique_classes; ++j) {
			printf("%5.2f ", v_m[j]);
		}
#endif
	}
	if ((n_perm > 0) && (per_perm != NULL)) {
		for (size_t j = 0; j < n_perm; ++j) {
			p_perm[j] = 100.0 * n_corr_perm[j] / (double) n;
		}
		qsort(p_perm, n_perm, sizeof(double), dcmpfunc);
		memcpy(per_perm, p_perm, n_perm * sizeof(double));
	}


#ifdef PRINT_DEBUG

	printf("\n\n");
	for(size_t i = 0; i < n_unique_classes; ++i)
	{
		printf("Item %zu: %zi\n", i, uclasses[i]);
	}
	printf("There are %zu unique items.\n", n_unique_classes);

	if (n_perm > 0) {
		size_t i_up = floor(n_perm * 0.975),
			i_low = floor(n_perm * 0.025),
			i_med = floor(n_perm * 0.5);
		printf("%zu %zu %zu\n", i_low, i_med, i_up);
		printf("Permutation correct: [%5.1f - %5.1f - %5.1f]\n", p_perm[i_low], p_perm[i_med], p_perm[i_up]);
	}
#endif

	if (knn_class) {
		for (size_t i = 0; i < n; ++i) {
			knn_class[i] = uclasses[_knn_class[i]];
		}
	}

	FREE_IF_NOT_NULL(uclasses);
	FREE_IF_NOT_NULL(nnclasses);
	FREE_IF_NOT_NULL(_knn_class);
	FREE_IF_NOT_NULL(v_m);
	FREE_IF_NOT_NULL(p_perm);
	FREE_IF_NOT_NULL(perm_ind);
	FREE_IF_NOT_NULL(n_corr_perm);
	FREE_IF_NOT_NULL(d_line);
	FREE_IF_NOT_NULL(weights);

	return 100.0 * num_correct / (double) n;
}

