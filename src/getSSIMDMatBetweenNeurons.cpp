/**
 * @file getSSIMDMatBetweenNeurons.cpp
 * \brief MATLAB wrapper for getSSIMDMatBetweenNeurons()
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
 *
 */

#include "mex.h"
#include "spikeTrainFun.hpp"

/**
 * MATLAB wrapper for getSSIMDMatBetweenNeurons().
 *
 * When compiled will produce a mex function which takes five arguments:
 * @param spikeTrains cell array of n vectors of doubles, representing spike trains from n neurons
 * @param evTimes vector of m doubles, representing times at which distances are taken
 * @param startWin beginning of window relative to evTimes (double scalar)
 * @param winLen window length (double scalar)
 * @param q Cost for moving a spike (temporal resolution; double scalar)
 * @return distance matrix, \f$m \times m \cdot n\f$, and cell array of m base spiketrains arrays, each n long
**/
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

	// check for proper number of arguments
	if ((nrhs < 5)) {
		mexErrMsgIdAndTxt("SSIMSToolbox:getSSIMDMatBetweenNeurons:nrhs", "This function expects 5 arguments.");
	}
	bool returnCell = false;
	if (nrhs > 5) {
		char cmd[64] = "asMatrix";
		if (mxGetString(prhs[5], cmd, sizeof(cmd))) {
			mexErrMsgIdAndTxt("SSIMSToolbox:getSSIMDMatBetweenNeurons:notAString", "Argument 6 (command) should be a string, and either be 'asMatrix' or 'asCell'.");
		}
		returnCell = (0 == strncmp(cmd, "asCell", 3));
	}

	if(nlhs < 1) {
		mexErrMsgIdAndTxt("SSIMSToolbox:getSSIMDMatBetweenNeurons:nlhs", "At least one output required: distance.");
	}
	// make sure the first input argument is type double
	const mxArray *spikeTrains = prhs[0];
	if( !mxIsCell(spikeTrains) || mxIsComplex(spikeTrains)) {
		mexErrMsgIdAndTxt("SSIMSToolbox:getSSIMDMatBetweenNeurons:notDouble","Expects cell array of spike trains as argument ar1.");
	}
	mwSize nSpikeTrains = mxGetNumberOfElements(spikeTrains);

	{
		mwIndex i;
		bool isNotDouble = false;
		for (i = 0; i < nSpikeTrains; ++i) {
			if (!mxIsDouble(mxGetCell(spikeTrains, i))) {
				isNotDouble = true;
				break;
			}
		}
		if (isNotDouble) {
			mexErrMsgIdAndTxt("SSIMSToolbox:getSSIMDMatBetweenNeurons:notDouble", "Element %i of spike train cell has class %s. We expect double!", i, mxGetClassName(mxGetCell(spikeTrains, i)));
		}
	}

	/* make sure the third input argument is scalar */
	if ( !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1])) {
		mexErrMsgIdAndTxt("SSIMSToolbox:getSSIMDMatBetweenNeurons:notDouble2", "Times (argument 2) must be real doubles.");
	}
	if (( !mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) || mxGetNumberOfElements(prhs[2])!=1)) {
		mexErrMsgIdAndTxt("SSIMSToolbox:getSSIMDMatBetweenNeurons:notScalar3","Window start must be scalar.");
	}
	if (( !mxIsDouble(prhs[3]) || mxIsComplex(prhs[3]) || mxGetNumberOfElements(prhs[3])!=1)) {
		mexErrMsgIdAndTxt("SSIMSToolbox:getSSIMDMatBetweenNeurons:notScalar4","Window length must be scalar.");
	}
	if (( !mxIsDouble(prhs[4]) || mxIsComplex(prhs[4]) || mxGetNumberOfElements(prhs[4])!=1)) {
		mexErrMsgIdAndTxt("SSIMSToolbox:getSSIMDMatBetweenNeurons:notScalar5","Cost q must be scalar.");
	}

	const mxArray *evTimes = prhs[1];
	double startWin = mxGetScalar(prhs[2]);
	double winLen = mxGetScalar(prhs[3]);
	double q = mxGetScalar(prhs[4]);

	mwSize nTimePoints = mxGetNumberOfElements(evTimes);
	mwSize nSpkBySpk = nSpikeTrains * (nSpikeTrains - 1) / 2;

	mxArray *baseSpk = NULL;
	if(nlhs > 1) {
		baseSpk = mxCreateCellMatrix(1, nTimePoints);
		plhs[1] = baseSpk;
	}

	SpkDataList_t tl(nSpikeTrains);

	mxArray * tmpSt;
	// Compile raw spike data
	for (SpkDataList_t::iterator s = tl.begin(); s != tl.end(); ++s) {
		tmpSt = mxGetCell(spikeTrains, s - tl.begin());
		*s = SpkData_t (mxGetPr(tmpSt), mxGetPr(tmpSt) + mxGetNumberOfElements(tmpSt));
	}
	sanitizeSpikeData(tl);

	if (returnCell) {
		Dist_t dMat(nSpkBySpk * nTimePoints);
		SpkDataList_t bspkl(nTimePoints * nSpikeTrains);
		getSSIMDMatBetweenNeurons(tl, SpkData_t(mxGetPr(evTimes), mxGetPr(evTimes) + mxGetNumberOfElements(evTimes)),
		startWin, winLen,  q,
		dMat.begin(), dMat.end(),
		bspkl.begin(), bspkl.end());

		/* Create output matrix */
		plhs[0] = mxCreateCellMatrix(nTimePoints, 1);
		mwIndex i_t = 0;
		for (Dist_t::iterator ii = dMat.begin(); ii != dMat.end(); ii += nSpkBySpk) {
			mxArray * nDst = mxCreateDoubleMatrix(nSpkBySpk, 1, mxREAL);
			std::copy(ii, ii + nSpkBySpk, mxGetPr(nDst));
			mxSetCell(plhs[0], i_t, nDst);
			++ i_t;
		}
		if (baseSpk) {
			for (mwIndex i_t = 0; i_t < nTimePoints; ++ i_t) {
				mxArray * nBspC = mxCreateCellMatrix(1, nSpikeTrains);
				for (mwIndex i_u = 0; i_u < nSpikeTrains; ++ i_u) {
					mwIndex ind = i_u + nSpikeTrains * i_t;
					mxArray * nBsp = mxCreateDoubleMatrix(1, bspkl[ind].size(), mxREAL);
					std::copy(bspkl[ind].begin(), bspkl[ind].end(), mxGetPr(nBsp));
					mxSetCell(nBspC, i_u, nBsp);
				}
				mxSetCell(baseSpk, i_t, nBspC);
			}
		}
	}
	else {
		Dist_t dMat(nSpkBySpk);
		SpkDataList_t bspkl(nSpikeTrains );
		double *evTimesDbl = mxGetPr(evTimes);
		/* Create output matrix */
		plhs[0] = mxCreateDoubleMatrix(nSpikeTrains, nSpikeTrains * nTimePoints, mxREAL);
		double *dmatDbl = mxGetPr(plhs[0]);
		for (mwIndex i_t = 0; i_t < nTimePoints; ++i_t) {
			getSSIMDMatBetweenNeurons(tl, SpkData_t(1, evTimesDbl[i_t]),
			startWin, winLen,  q,
			dMat.begin(), dMat.end(),
			bspkl.begin(), bspkl.end());

			copyPDistToSquare(dMat, dmatDbl + nSpikeTrains * nSpikeTrains * i_t);
			if (baseSpk) {
				mxArray * nBspC = mxCreateCellMatrix(1, nSpikeTrains);
				for (mwIndex i_u = 0; i_u < nSpikeTrains; ++ i_u) {
					mxArray * nBsp = mxCreateDoubleMatrix(1, bspkl[i_u].size(), mxREAL);
					std::copy(bspkl[i_u].begin(), bspkl[i_u].end(), mxGetPr(nBsp));
					mxSetCell(nBspC, i_u, nBsp);
				}
				mxSetCell(baseSpk, i_t, nBspC);
			}

		}
	}
}
