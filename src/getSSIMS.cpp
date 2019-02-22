/**
 * @file getSSIMS.cpp
 * @brief MATLAB wrapper for TSNERunner class, to run tSNE on some data, (mostly) replicating getSSIMS.m
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
#include "armaMex.hpp"
#include "spikeTrainFun.hpp"
#include "tSNERunnerD.hpp"
#include "tSNERunnerV.hpp"

/**	MATLAB wrapper for TSNERunner, replicating getSSIMS.m.
 *
 * 	When compiled will produce a mex function which takes three arguments:
 * 	@param spikeTrains cell array of n vectors of doubles, representing spike trains from n neurons
 * 	@param evTimes vector of m doubles, representing times at which distances are taken
 * 	@param winLen window length (double scalar, default = 1.0s), covering spikes in window (evTimes, evTimes + winLen)
 * 	@param q Cost for moving a spike (temporal resolution; double scalar, default = 10.0)
 * 	@param outDimensions number of dimensions the dimensionality-reduced data should have. This number should be less than
 * 		number of columns of baseData (optional, default = 3)
 * 	@param perplexity Sets perplexity parameter for tSNE (optional, default = 30.0)
 * 	@param direct_tSNE: if true, and only one spike train given, then tSNE will be performed on the distance matrix.
 * 		If false, PCA will be performed first (optional, default = false)
 *
 * 	@return 	tSNE-transformed data,
 * 				tSNE-transformation matrix to transform more data,
 * 				baseSpikeTrains,
 * 				baseDMat,
 * 				kld Kullback-Leibler-divergence between initial and low-dimensional probability distributions
**/
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

	// check for proper number of arguments
	if (nrhs < 2) {
		mexErrMsgIdAndTxt("SSIMSToolbox:runtSNE:nrhs", "This function expects at least two arguments.");
	}

	// make sure the first input argument is type double
	const mxArray *spikeTrains = prhs[0];

	if( !mxIsCell(spikeTrains) || mxIsComplex(spikeTrains)) {
		mexErrMsgIdAndTxt("SSIMSToolbox:getSSIMS:notDouble", "Expects cell array of spike trains as argument 1.");
	}
	mwSize nSpikeTrains = mxGetNumberOfElements(spikeTrains);

	{
		bool isNotDouble = false;
		mwIndex i;
		for (i = 0; i < nSpikeTrains; ++i) {
			if (!mxIsDouble(mxGetCell(prhs[0], i))) {
				isNotDouble = true;
				break;
			}
		}
		if (isNotDouble) {
			mexErrMsgIdAndTxt("SSIMSToolbox:getSSIMS:notDouble","Element %i of spike train cell has class %s. We expect double!", i, mxGetClassName(mxGetCell(prhs[0], i)));
		}
	}

	// make sure other arguments have correct format
	const mxArray *evTimes = prhs[1];
	if ( !mxIsDouble(evTimes) || mxIsComplex(evTimes)) {
		mexErrMsgIdAndTxt("SSIMSToolbox:getSSIMS:notDouble2", "Times (argument 2) must be real doubles.");
	}
	mwSize nTimePoints = mxGetNumberOfElements(evTimes);


	double winLen = 1.0;
	if (nrhs > 2) {
		if ( !mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) || mxGetNumberOfElements(prhs[2]) != 1) {
			mexErrMsgIdAndTxt("SSIMSToolbox:getSSIMS:notScalar3", "Window length must be scalar.");
		}
		winLen = mxGetScalar(prhs[2]);
	}

	double q = 10.0;
	if (nrhs > 3) {
		if ( !mxIsDouble(prhs[3]) || mxIsComplex(prhs[3]) || mxGetNumberOfElements(prhs[3]) !=1 ) {
			mexErrMsgIdAndTxt("SSIMSToolbox:getSSIMS:notScalar4", "Cost q must be scalar.");
		}
		q = mxGetScalar(prhs[3]);
	}

	mwSize nDims = 3;
	if (nrhs > 4) {
		if ( !mxIsDouble(prhs[4]) || mxIsComplex(prhs[4]) || mxGetNumberOfElements(prhs[4]) !=1 ) {
			mexErrMsgIdAndTxt("SSIMSToolbox:getSSIMS:notScalar5", "Out dimensions must be scalar.");
		}
		nDims = mxGetScalar(prhs[4]);

	}

	double perplexity = 30.0;
	if (nrhs > 5) {
		if ( !mxIsDouble(prhs[5]) || mxIsComplex(prhs[5]) || mxGetNumberOfElements(prhs[5]) != 1) {
			mexErrMsgIdAndTxt("SSIMSToolbox:getSSIMS:notScalar6", "Perplexity must be scalar.");
		}
		perplexity = mxGetScalar(prhs[5]);
	}

	bool direct_tSNE = false;
	if (nrhs > 6) {
		if ( (!mxIsDouble(prhs[6]) && !mxIsClass(prhs[6], "logical")) || mxIsComplex(prhs[6]) || mxGetNumberOfElements(prhs[6]) != 1) {
			mexErrMsgIdAndTxt("SSIMSToolbox:getSSIMS:notScalar7", "Direct_tSNE must be scalar.");
		}
		direct_tSNE = mxGetScalar(prhs[6]) != 0;
	}

	SpkDataList_t baseSpkTr;
	mxArray *baseSpk = NULL; // output base spike trains

	if(nlhs > 2) {
		baseSpk = mxCreateCellMatrix(1, nSpikeTrains);
		plhs[2] = baseSpk;
		baseSpkTr.resize(nSpikeTrains * nTimePoints);
	}

	SpkDataList_t tl(nSpikeTrains);
	mxArray * tmpSt;
	// Compile raw spike data
	for (SpkDataList_t::iterator s = tl.begin(); s != tl.end(); ++s) {
		tmpSt = mxGetCell(spikeTrains, s - tl.begin());
		*s = SpkData_t (mxGetPr(tmpSt), mxGetPr(tmpSt) + mxGetNumberOfElements(tmpSt));
	}
	sanitizeSpikeData(tl);


	arma::mat aBaseData = getSSIMDMatBetweenTimePoints( tl, SpkData_t(mxGetPr(evTimes), mxGetPr(evTimes) + mxGetNumberOfElements(evTimes)),
	0.0, winLen, q, baseSpkTr.begin(), baseSpkTr.end());

	TSNERunner *tSNE;

	if (direct_tSNE && (nSpikeTrains == 1)) {
		tSNE = new TSNERunnerD(aBaseData, nDims, perplexity);
	}
	else {
		tSNE = new TSNERunnerV(aBaseData, nDims, perplexity);
	}

	double kld = tSNE->run();

	if (nlhs > 0) {
		const arma::mat ydata = tSNE->getTSNEData();
		plhs[0] = armaCreateMxMatrix(ydata.n_rows, ydata.n_cols);
		armaSetPr(plhs[0], ydata);
	}
	if (nlhs > 1) {
		if (direct_tSNE) {
			plhs[1] = mxCreateDoubleMatrix(0, 0, mxREAL);
		}
		else {
			TSNERunnerV * tSNEV = static_cast<TSNERunnerV*>(tSNE);
			const arma::mat transform = tSNEV->getTSNETransform();
			plhs[1] = armaCreateMxMatrix(transform.n_rows, transform.n_cols);
			armaSetPr(plhs[1], transform);
		}
	}
	if (nlhs > 2) {
		for (mwIndex i = 0; i < nSpikeTrains; ++i) {
			mxArray *tmpBaseCell = mxCreateCellMatrix(1, nTimePoints);
			for (mwIndex j = 0; j < nTimePoints; ++j) {
				mwIndex k = i * nTimePoints + j;
				mxArray * nBsp = mxCreateDoubleMatrix(1, baseSpkTr[k].size(), mxREAL);
				std::copy(baseSpkTr[k].begin(), baseSpkTr[k].end(), mxGetPr(nBsp));
				mxSetCell(tmpBaseCell, j, nBsp);
			}
			mxSetCell(baseSpk, i, tmpBaseCell);
		}
	}
	if (nlhs > 3) {
		plhs[3] = armaCreateMxMatrix(aBaseData.n_rows, aBaseData.n_cols);
		armaSetPr(plhs[3], aBaseData);
	}
	if (nlhs > 4) {
		plhs[4] = mxCreateDoubleScalar(kld);
	}
	delete tSNE;
}
