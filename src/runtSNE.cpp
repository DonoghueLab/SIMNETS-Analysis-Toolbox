/**
 * @file runtSNE.cpp
 * @brief MATLAB wrapper for TSNERunner class, to run tSNE on some data
 *
 * @author  Jonas Zimmermann <jonas_zimmermann@brown.edu>
 * @version 0.0.1
 *
 * @section LICENSE
 *
 * This file is part of SSIMS Toolbox.
 *
 * SSIMS Toolbox is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SSIMS Toolbox is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with SSIMS Toolbox.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include "mex.h"
#include "tSNERunnerV.hpp"
#include "armaMex.hpp"

/**	MATLAB wrapper for TSNERunner.
 *
 * 	When compiled will produce a mex function which takes three arguments:
 * 	@param baseData matrix containing data to be transformed. Rows contain observations, columns contain features
 * 	@param outDimensions number of dimensions the dimensionality-reduced data should have. This number should be less than
 * 		number of columns of baseData (optional, default = 2)
 * 	@param perplexity Sets perplexity parameter for tSNE (optional, default = 30)
 *
 * 	@return tSNE-transformed data, as well as a tSNE-transformation matrix to transform more data, and also the
 * 		Kullback-Leibler-divergence between initial and low-dimensional probability distributions
**/
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    try {
        // check for proper number of arguments
        if (nrhs < 1) {
            mexErrMsgIdAndTxt("SSIMSToolbox:runtSNE:nrhs", "This function expects at least one argument.");
        }
        if (mxIsComplex(prhs[0]) || mxIsCell(prhs[0]) || !mxIsDouble(prhs[0])) {
            mexErrMsgIdAndTxt("SSIMSToolbox:runtSNE:nrhs", "First argument is expected to be a real double matrix.");
        }
        const mxArray *baseData = prhs[0];

        // Set number of output dimensions
        mwSize nDims = 2;
        if (nrhs > 1) {
            nDims = (mwSize) mxGetScalar(prhs[1]);
        }

        // Set perplexity
        double perplexity = 30.0;
        if (nrhs > 2) {
            perplexity = mxGetScalar(prhs[2]);
        }

        if(nlhs < 1) {
            mexErrMsgIdAndTxt("SSIMSToolbox:runtSNE:nlhs", "At least one output required: tSNE output.");
        }

        arma::mat aBaseData = armaGetPr(baseData);
        TSNERunnerV tSNE(aBaseData, nDims, perplexity);
        double kld = tSNE.run();

        if (nlhs > 0) {
            const arma::mat ydata = tSNE.getTSNEData();
            plhs[0] = armaCreateMxMatrix(ydata.n_rows, ydata.n_cols);
            armaSetPr(plhs[0], ydata);
        }
        if (nlhs > 1) {
            const arma::mat transform = tSNE.getTSNETransform();
            plhs[1] = armaCreateMxMatrix(transform.n_rows, transform.n_cols);
            armaSetPr(plhs[1], transform);
        }
        if (nlhs > 2) {
            plhs[2] = mxCreateDoubleScalar(kld);
        }
    }
    catch(std::exception e){
        mexErrMsgIdAndTxt("SSIMS:exception", "Running tSNE: %s", e.what());
    }
}
