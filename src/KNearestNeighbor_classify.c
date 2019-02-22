#include "mex.h"
#include <math.h>
#include <string.h>
#include <time.h>
#include <stdint.h>
#include "knnHelpers.h"

// double doKnnClassification(const double *dist_mat, const size_t n, const long long *true_class, const size_t k, const size_t n_perm, long long *knn_class, double *vote_matrix, double *per_perm);


void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
	// input to KNN classifier
	size_t n = 0;					// number of points
	double *data = NULL;			// raw data
	size_t n_cols = 0;				// number of columns in data matrix
	size_t n_rows = 0;				// number of rows in data matrix
	long long *id = NULL;			// group id vector (size n!)
	size_t k = 1;					// k nearest neighbours will be looked at
	size_t n_perm = 0;				// number of permutations for test

	// derived from id
	ssize_t n_unique_ids = 0;

	double *dist_mat_tril = NULL;		// distance matrix between points (n * n)
	// output from KNN classifier
	long long *knn_id = NULL;		// will hold ID result from KNN classifier
	double *vote_matrix = NULL;		// matrix that holds n * (number of unique IDs) how IDs voted
	double *per_perm = NULL;		// percent correct (size n_perm) from permutation test
	double percent_correct = 0.0;	// percent correct from KNN classification with cross-validation

	/* check for proper number of arguments */
	if ((nrhs < 3) || (nrhs > 4)) {
		mexErrMsgIdAndTxt("SSIMS:KNearestNeighbor_classify:nrhs","This function expects 3 or 4 inputs.");
	}
	/* make sure the first input argument is type double */
	if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0])) {
		mexErrMsgIdAndTxt("SSIMS:KNearestNeighbor_classify:notDouble1","Expects real double distance matrix as first argument.");
	}
	/* make sure the second input argument is type double */
	if ((!mxIsNumeric(prhs[1]) || mxIsComplex(prhs[1])) ) {
		mexPrintf("arg 2 is of type: %s\n", mxGetClassName(prhs[1]));
		mexErrMsgIdAndTxt("SSIMS:KNearestNeighbor_classify:notDouble2","Expects array of IDs be an integer or real double type.");
	}

	if ( (!mxIsNumeric(prhs[2])) || mxIsComplex(prhs[2]) || (mxGetNumberOfElements(prhs[2]) != 1 ) ) {
		mexErrMsgIdAndTxt("SSIMS:KNearestNeighbor_classify:notDouble3","Expects k to be a numeric scalar.");
	}

	if ((nrhs > 3) && ( (!mxIsNumeric(prhs[3])) || mxIsComplex(prhs[3]) || (mxGetNumberOfElements(prhs[3]) != 1 ) )) {
		mexErrMsgIdAndTxt("SSIMS:KNearestNeighbor_classify:notDouble4","Expects n_perm to be a numeric scalar.");
	}

	k = (size_t) mxGetScalar(prhs[2]);
	if (nrhs > 3) {
		n_perm = (size_t) mxGetScalar(prhs[3]);
	}
	n = mxGetNumberOfElements(prhs[1]);

	if ( k > n - 1 ) {
		mexWarnMsgIdAndTxt("SSIMS:KNearestNeighbor_classify:kGreaterThanN","Expects k (given: %zu) to be strictly less than n (%zu). We will set k = %zu.", k, n, n - 1 );
		k = n - 1;
	}

	if ((mxGetM(prhs[0]) != n ) && (mxGetN(prhs[0]) != n )) {
		mexErrMsgIdAndTxt("SSIMS:KNearestNeighbor_classify:dataMalformed","Neither row nor column number of data equals number of IDs (%zu).", n);
	}
	n_cols = mxGetN(prhs[0]);
	n_rows = mxGetM(prhs[0]);

	data = mxGetPr(prhs[0]);
	dist_mat_tril = mxMalloc(n * (n - 1) / 2 * sizeof(double));
	if ( n == n_cols ) {
		Rpdist_c(data, n_rows, n_cols, dist_mat_tril);
	} else {
		mexWarnMsgIdAndTxt("SSIMS:KNearestNeighbor_classify:dataTransposed","Data has n = %zu observations, but matrix has %zu columns. We will therefore transpose it. This may not be what you want!", n, n_cols );
		Rpdist(data, n_rows, n_cols, dist_mat_tril);
	}


	// initialize random number generator
	srand((unsigned) time(NULL));

	mxClassID id_class = mxGetClassID(prhs[1]);
	id = mxMalloc(n * sizeof(long long));
	if ((id_class == mxINT64_CLASS) || (id_class == mxUINT64_CLASS))
		memcpy(id, mxGetData(prhs[1]), n * sizeof(long long));
	else {
		void *data = mxGetData(prhs[1]);
		switch (id_class)
		{
			case mxINT8_CLASS:
			case mxUINT8_CLASS:
				for (size_t i = 0; i < n; i++) {
					id[i] = (long long) ((int8_t *) data)[i];
				}
				break;
			case mxINT16_CLASS:
			case mxUINT16_CLASS:
				for (size_t i = 0; i < n; i++) {
					id[i] = (long long)((int16_t *) data)[i];
				}
				break;
			case mxINT32_CLASS:
			case mxUINT32_CLASS:
				for (size_t i = 0; i < n; i++) {
					id[i] = (long long)((int *) data)[i];
				}
				break;
			case mxSINGLE_CLASS:
				if (sizeof(float) != 4) {
					mexErrMsgIdAndTxt("SSIMS:KNearestNeighbor_classify:weirdFloat","Oh no. float is not a 32bit type. That means automatic handling will break. Would you mind passing group IDs as integers, e.g. int64(id)?");
				}
				for (size_t i = 0; i < n; i++) {
					id[i] = (long long)((int32_t *) data)[i];
				}
				break;
			case mxDOUBLE_CLASS:
				if (sizeof(double) == sizeof(long long)) {
					for (size_t i = 0; i < n; i++) {
						id[i] = (long long)((long long*) data)[i];
					}
				} else if (sizeof(double) == sizeof(int)) {
					for (size_t i = 0; i < n; i++) {
						id[i] = (long long)((int *) data)[i];
					}
				} else {
					mexErrMsgIdAndTxt("SSIMS:KNearestNeighbor_classify:weirdDouble","Oh no. double is neither the size of int nor long long. That means automatic handling will break. Would you mind passing group IDs as integers, e.g. int64(id)?");
				}
				break;
			default:
				mexErrMsgIdAndTxt("SSIMS:KNearestNeighbor_classify:notNumeric","id must be (u)int(8|16|32|64) or single/double.");
				break;
		}
	}

	plhs[0] = mxCreateNumericMatrix(n, 1, id_class, mxREAL);

	if (nlhs > 2) {
		long long *ulist = mxMalloc(n * sizeof(long long));
		n_unique_ids = findUniqueItemsLL(id, n, ulist);
		plhs[2] = mxCreateDoubleMatrix(n_unique_ids, n, mxREAL);
		vote_matrix = mxGetPr(plhs[2]);
		mxFree(ulist);
		ulist = NULL;
	}
	if (nlhs > 3) {
		plhs[3] = mxCreateDoubleMatrix(n_perm, 1, mxREAL);
		per_perm = mxGetPr(plhs[3]);
	}

	knn_id = mxMalloc(n * sizeof(long long));
	percent_correct = doKnnClassification(dist_mat_tril, n, id, k, n_perm, knn_id, vote_matrix, per_perm);

	if (nlhs > 1) {
		plhs[1] = mxCreateDoubleScalar(percent_correct);
	}

	if ((id_class == mxINT64_CLASS) || (id_class == mxUINT64_CLASS))
		memcpy(knn_id, mxGetData(plhs[0]), n * sizeof(long long));
	else {
		void *data = mxGetData(plhs[0]);
		switch (id_class)
		{
			case mxINT8_CLASS:
			case mxUINT8_CLASS:
				for (size_t i = 0; i < n; i++) {
					((int8_t *)data)[i] = knn_id[i];
				}
				break;
			case mxINT16_CLASS:
			case mxUINT16_CLASS:
				for (size_t i = 0; i < n; i++) {
					((int16_t *) data)[i] = knn_id[i];
				}
				break;
			case mxINT32_CLASS:
			case mxUINT32_CLASS:
				for (size_t i = 0; i < n; i++) {
					((int *) data)[i] = knn_id[i];
				}
				break;
			case mxSINGLE_CLASS:
				for (size_t i = 0; i < n; i++) {
					((float *) data)[i] = *(float *) &(knn_id[i]);
				}
				break;
			case mxDOUBLE_CLASS:
				for (size_t i = 0; i < n; i++) {
					((double *)data)[i] = *(double *) &(knn_id[i]);
				}
				break;
			default:
				mexErrMsgIdAndTxt("SSIMS:KNearestNeighbor_classify:notNumeric","id must be (u)int(8|16|32|64) or single/double.");
				break;
		}
	}


	// clean up
	mxFree(dist_mat_tril);
	mxFree(id);
	mxFree(knn_id);
}
