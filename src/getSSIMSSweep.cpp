/**
 * @file getSSIMSSweep.cpp
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
#include "ssimsHelper.hpp"

/**	MATLAB wrapper for doSsimsSweep.
 *
 * 	When compiled will produce a mex function which takes eight arguments:
 * 	@param spikeTrains cell array of n vectors of doubles, representing spike trains from n neurons
 * 	@param evTimes vector of m doubles, representing times at which distances are taken
 *	@param evOffsets vector of n_o doubles, representing offset times for events over which SSIMS will be looped
 * 	@param winLen window length (double scalar, default = 1.0s), covering spikes in window (evTimes, evTimes + winLen) (or vector of such times)
 * 	@param q Vector of costs for moving a spike (temporal resolution; double scalar, default = 10.0), to be looped over
 * 	@param outDimensions number of dimensions the dimensionality-reduced data should have. This number should be less than
 * 		number of columns of baseData (optional, default = 3), or vector; will be looped over
 * 	@param perplexity Sets perplexity parameter for tSNE (optional, default = 30.0), or vector; will be looped over
 * 	@param direct_tSNE: if true, and only one spike train given, then tSNE will be performed on the distance matrix.
 * 		If false, PCA will be performed first (optional, default = false)
 *
 * 	@return 	tSNE-transformed data, cell of size |evOffsets| * |winLen| * |q| * |outDimensions| * |perplexity|
**/
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

	// check for proper number of arguments
	if (nrhs < 2) {
		mexErrMsgIdAndTxt("SSIMSToolbox:runtSNE:nrhs", "This function expects at least two arguments.");
	}

	// make sure the first input argument is a cell of type double
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

	double *evOffsets, evOffset = 0;
	mwSize nEventOffsets = 1;
	if (nrhs > 2) {
		evOffsets = mxGetPr(prhs[2]);
		nEventOffsets = mxGetNumberOfElements(prhs[2]);
		if ( !mxIsDouble(prhs[2]) || mxIsComplex(prhs[2])) {
			mexErrMsgIdAndTxt("SSIMSToolbox:getSSIMS:notDouble2", "Event Offsets (argument 3) must be real doubles.");
		}
	}
	else {
		evOffsets = &evOffset;
	}

	double *winLens, winLen = 1.0;
	mwSize nWinLens = 1;

	if (nrhs > 3) {
		if ( !mxIsDouble(prhs[3]) || mxIsComplex(prhs[3])) {
			mexErrMsgIdAndTxt("SSIMSToolbox:getSSIMS:notScalar3", "Window lengths must be doubles.");
		}
		winLens = mxGetPr(prhs[3]);
		nWinLens = mxGetNumberOfElements(prhs[3]);
	}
	else {
		winLens = &winLen;
	}

	double *qs, q = 10.0;
	mwSize nQs = 1;
	if (nrhs > 4) {
		if ( !mxIsDouble(prhs[4]) || mxIsComplex(prhs[4])) {
			mexErrMsgIdAndTxt("SSIMSToolbox:getSSIMS:notScalar4", "Costs q must be doubles.");
		}
		qs = mxGetPr(prhs[4]);
		nQs = mxGetNumberOfElements(prhs[4]);
	}
	else {
		qs = &q;
	}

	double * nDims, nDim = 3;
	mwSize nNDims = 1;
	if (nrhs > 5) {
		if ( !mxIsDouble(prhs[5]) || mxIsComplex(prhs[5])) {
			mexErrMsgIdAndTxt("SSIMSToolbox:getSSIMS:notScalar5", "Out dimensions must be a vector of integers.");
		}
		nDims = mxGetPr(prhs[5]);
		nNDims = mxGetNumberOfElements(prhs[5]);

	}
	else {
		nDims = &nDim;
	}

	double * perplexities, perplexity = 30.0;
	mwSize nPerplexities = 1;
	if (nrhs > 6) {
		if ( !mxIsDouble(prhs[6]) || mxIsComplex(prhs[6])) {
			mexErrMsgIdAndTxt("SSIMSToolbox:getSSIMS:notScalar6", "Perplexities must be doubles.");
		}
		perplexities = mxGetPr(prhs[6]);
		nPerplexities = mxGetNumberOfElements(prhs[6]);
	}
	else {
		perplexities = &perplexity;
	}


	bool direct_tSNE = false;
	if (nrhs > 7) {
		if ( (!mxIsDouble(prhs[7]) && !mxIsClass(prhs[7], "logical")) || mxIsComplex(prhs[7]) || mxGetNumberOfElements(prhs[7]) != 1) {
			mexErrMsgIdAndTxt("SSIMSToolbox:getSSIMS:notScalar7", "Direct_tSNE must be scalar.");
		}
		direct_tSNE = mxGetScalar(prhs[7]) != 0;
	}

	SpkDataList_t baseSpkTr;

	SpkDataList_t tl(nSpikeTrains);
	mxArray * tmpSt;
	// Compile raw spike data
	for (SpkDataList_t::iterator s = tl.begin(); s != tl.end(); ++s) {
		tmpSt = mxGetCell(spikeTrains, s - tl.begin());
		*s = SpkData_t (mxGetPr(tmpSt), mxGetPr(tmpSt) + mxGetNumberOfElements(tmpSt));
	}
	sanitizeSpikeData(tl);

	if (nlhs > 0) {
		direct_tSNE = direct_tSNE && (nSpikeTrains == 1);
		// Run SSIMS sweep
		std::vector<arma::mat> outSSIMS = doSsimsSweep(tl, SpkData_t(mxGetPr(evTimes), mxGetPr(evTimes) + mxGetNumberOfElements(evTimes)), SpkData_t(evOffsets, evOffsets + nEventOffsets), SpkData_t(winLens, winLens + nWinLens), SpkData_t(qs, qs + nQs), std::vector<size_t>(nDims, nDims + nNDims), std::vector<double>(perplexities, perplexities + nPerplexities), direct_tSNE);

		mwSize resSize[] = {nEventOffsets, nWinLens, nQs, nNDims, nPerplexities};
		plhs[0] = mxCreateCellArray(5, resSize);

		// Prepare output
		for (mwSize i_p = 0; i_p < nPerplexities; ++i_p) {
			for (mwSize i_d = 0; i_d < nNDims; ++i_d) {
				for (mwSize i_q = 0; i_q < nQs; ++i_q) {
					for (mwSize i_wl = 0; i_wl < nWinLens; ++i_wl) {
						for (mwSize i_eo = 0; i_eo < nEventOffsets; ++i_eo) {
							mwIndex indA[] = {i_eo, i_wl, i_q, i_d, i_p};
							mwIndex ind = mxCalcSingleSubscript(plhs[0], 5, indA);
							mxArray * SSIMSr = armaCreateMxMatrix(outSSIMS[ind].n_rows, outSSIMS[ind].n_cols);
							armaSetPr(SSIMSr, outSSIMS[ind]);

 							mxSetCell(plhs[0], ind, SSIMSr);
						}
					}
				}
			}
		}
	}
}
