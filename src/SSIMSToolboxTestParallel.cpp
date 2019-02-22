#include "mex.h"
#if defined(_OPENMP)
#include <omp.h>
#endif
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
#if defined(_OPENMP)
	if (nrhs > 0) {
		omp_set_num_threads(mxGetScalar(prhs[0]));
	}
	{
        mexPrintf("If parallel processing is enabled, you should see the line\n\"Hello world from thread x\", where 'x' represents a thread number.\nThere will be %i threads in parallel.\n\n", omp_get_max_threads());
		#pragma omp parallel
		{
			mexPrintf("Hello world from thread %i.\n", omp_get_thread_num());
		}
	}
#else
	mexPrintf("SSIMS Toolbox was not compiled with support for parallelization.\n\n");
#endif
}
