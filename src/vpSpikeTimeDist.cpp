/**
 * @file vpSpikeTimeDist.cpp
 * \brief MATLAB wrapper to get VP spike train distances
 *
 * @author  Jonas Zimmermann <jonas_zimmermann@brown.edu>
 * @version 0.5
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
**/

#include "mex.h"
#include "spikeTrainFun.hpp"

/**
 * MATLAB wrapper to get VP spike train distances.
 *
 * When compiled will produce a mex function which takes three arguments:
 * @param spkT1 is a single spike train (vector of spike times [double])
 * @param spkT2 is a single spike train (vector of spike times [double])
 * @param q is a double scalar, the cost of moving a spike
 * @return spike train distance, a double scalar
**/
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	// check for proper number of arguments
	if (nrhs < 3) {
		mexErrMsgIdAndTxt("SSIMSToolbox:vpSpikeTimeDist:nrhs", "This function expects at least 3 inputs.");
	}
	SpkData_t spkd1, spkd2;
	SpikeTrain spkt1, spkt2;

	if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) || !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1])) {
		mexErrMsgIdAndTxt("SSIMSToolbox:vpSpikeTimeDist:spkNotDouble", "Spike times must be real double!");
	}
	spkd1 = SpkData_t(mxGetPr(prhs[0]), mxGetPr(prhs[0]) + mxGetNumberOfElements(prhs[0]));
	spkd2 = SpkData_t(mxGetPr(prhs[1]), mxGetPr(prhs[1]) + mxGetNumberOfElements(prhs[1]));
	sanitizeSpikeData(spkd1);
	sanitizeSpikeData(spkd2);
	spkt1 = SpikeTrain(spkd1);
	spkt2 = SpikeTrain(spkd2);

	if (!mxIsDouble(prhs[2]) || mxIsComplex(prhs[2])) {
		mexErrMsgIdAndTxt("SSIMSToolbox:vpSpikeTimeDist:qNotDouble", "Spike times must be real double!");
	}
	double q = mxGetScalar(prhs[2]);



	mxArray *dist = NULL;
	if(nlhs > 0) {
		dist = mxCreateDoubleScalar(spkt1.distVPTo(spkt2, q));
		plhs[0] = dist;

	}

}
