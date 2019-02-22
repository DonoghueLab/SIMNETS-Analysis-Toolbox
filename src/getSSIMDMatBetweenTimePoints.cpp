/**
 * @file getSSIMDMatBetweenTimePoints.cpp
 * \brief MATLAB wrapper for getSSIMDMatBetweenTimePoints()
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
#include <cstring>

#include "spikeTrainFun.hpp"

/**
 * MATLAB wrapper for getSSIMDMatBetweenTimePoints().
 *
 * When compiled will produce a mex function which takes five arguments:
 * @param spikeTrains cell array of n vectors of doubles, representing spike trains from n neurons
 * @param evTimes vector of m doubles, representing times at which distances are taken
 * @param winLen window length (double scalar)
 * @param q Cost for moving a spike (temporal resolution; double scalar)
 * @param cmd optional string argument to determine whether the distances will be returned as
 *  double matrix ('asMatrix'; default) or as a cell array ('asCell')
 * @return distance matrix, \f$n \times n \cdot m\f$, and cell array of n base spiketrains arrays, each m long
**/
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	mwSize nTpByTp;              // and their sies

	// check for proper number of arguments
	if (nrhs < 4) {
		mexErrMsgIdAndTxt("SSIMSToolbox:getSSIMDMatBetweenTimePoints:nrhs", "This function expects 4 inputs.");
	}

	bool returnCell = false;
	if (nrhs > 4) {
		char cmd[64] = "asMatrix";
		if (mxGetString(prhs[4], cmd, sizeof(cmd))) {
			mexErrMsgIdAndTxt("SSIMSToolbox:getSSIMDMatBetweenTimePoints:notAString", "Argument 5 (command) should be a string, and either be 'asMatrix' or 'asCell'.");
		}
		returnCell = (0 == strncmp(cmd, "asCell", 3));
	}

	if(nlhs < 1) {
		mexErrMsgIdAndTxt("SSIMSToolbox:getSSIMDMatBetweenTimePoints:nlhs", "At least one output required: distance.");
	}
	// make sure the first input argument is type double
	const mxArray *spikeTrains = prhs[0];

	if( !mxIsCell(spikeTrains) || mxIsComplex(spikeTrains)) {
		mexErrMsgIdAndTxt("SSIMSToolbox:getSSIMDMatBetweenTimePoints:notDouble", "Expects cell array of spike trains as argument 1.");
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
			mexErrMsgIdAndTxt("SSIMSToolbox:getSSIMDMatBetweenTimePoints:notDouble","Element %i of spike train cell has class %s. We expect double!", i, mxGetClassName(mxGetCell(spikeTrains, i)));
		}
	}
	// make sure other arguments have correct format
	const mxArray *evTimes = prhs[1];
	if ( !mxIsDouble(evTimes) || mxIsComplex(evTimes)) {
		mexErrMsgIdAndTxt("SSIMSToolbox:getSSIMDMatBetweenTimePoints:notDouble2", "Times (argument 2) must be real doubles.");
	}
	mwSize nTimePoints = mxGetNumberOfElements(evTimes);

	double startWin = 0.0;

	if ( !mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) || mxGetNumberOfElements(prhs[2]) != 1) {
		mexErrMsgIdAndTxt("SSIMSToolbox:getSSIMDMatBetweenTimePoints:notScalar3","Window length must be scalar.");
	}
	double winLen = mxGetScalar(prhs[2]);

	if ( !mxIsDouble(prhs[3]) || mxIsComplex(prhs[3]) || mxGetNumberOfElements(prhs[3]) != 1) {
		mexErrMsgIdAndTxt("SSIMSToolbox:getSSIMDMatBetweenTimePoints:notScalar4","Cost q must be scalar.");
	}
	double q = mxGetScalar(prhs[3]);

	nTpByTp = nTimePoints * (nTimePoints - 1) / 2;

	mxArray *baseSpk = NULL;
	if(nlhs > 1) {
		baseSpk = mxCreateCellMatrix(nSpikeTrains, 1);
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
		Dist_t dMat(nTpByTp * nSpikeTrains);
		SpkDataList_t bspkl(nTimePoints * nSpikeTrains);

		getSSIMDMatBetweenTimePoints(tl, SpkData_t(mxGetPr(evTimes), mxGetPr(evTimes) + mxGetNumberOfElements(evTimes)),
		startWin, winLen,  q,
		dMat.begin(), dMat.end(),
		bspkl.begin(), bspkl.end());
		/* Create output matrix */
		plhs[0] = mxCreateCellMatrix(nSpikeTrains, 1);
		mwIndex i_t = 0;
		for (Dist_t::iterator ii = dMat.begin(); ii != dMat.end(); ii += nTpByTp) {
			mxArray * nDst = mxCreateDoubleMatrix(nTpByTp, 1, mxREAL);
			std::copy(ii, ii + nTpByTp, mxGetPr(nDst));
			mxSetCell(plhs[0], i_t, nDst);
			++ i_t;
		}

		mwIndex ind;
		if (baseSpk) {
			for (mwIndex i_u = 0; i_u < nSpikeTrains; ++i_u) {
				mxArray * nBspC = mxCreateCellMatrix(nTimePoints, 1);
				for (mwIndex r = 0; r < nTimePoints; ++r) {
					ind = r + nTimePoints * i_u;
					mxArray * nBsp = mxCreateDoubleMatrix(1, bspkl[ind].size(), mxREAL);
					std::copy(bspkl[ind].begin(), bspkl[ind].end(), mxGetPr(nBsp));
					mxSetCell(nBspC, r, nBsp);
				}
				mxSetCell(baseSpk, i_u, nBspC);
			}
		}
	}
	else {
		Dist_t dMat(nTpByTp);
		SpkDataList_t bspkl(nTimePoints );
		/* Create output matrix */
		plhs[0] = mxCreateDoubleMatrix(nTimePoints, nTimePoints * nSpikeTrains, mxREAL);
		double *dmatDbl = mxGetPr(plhs[0]);
		for (mwIndex i_u = 0; i_u < nSpikeTrains; ++i_u) {
			getSSIMDMatBetweenTimePoints(SpkDataList_t(&(tl[i_u]), (&(tl[i_u]))+1), SpkData_t(mxGetPr(evTimes), mxGetPr(evTimes) + mxGetNumberOfElements(evTimes)),
			startWin, winLen,  q,
			dMat.begin(), dMat.end(),
			bspkl.begin(), bspkl.end());

			copyPDistToSquare(dMat, dmatDbl + nTimePoints * nTimePoints * i_u);

			if (baseSpk) {
				mxArray * nBspC = mxCreateCellMatrix(nTimePoints, 1);
				for (mwIndex i_t = 0; i_t < nTimePoints; ++ i_t) {
					mxArray * nBsp = mxCreateDoubleMatrix(1, bspkl[i_t].size(), mxREAL);
					std::copy(bspkl[i_t].begin(), bspkl[i_t].end(), mxGetPr(nBsp));
					mxSetCell(nBspC, i_t, nBsp);
				}
				mxSetCell(baseSpk, i_u, nBspC);
			}
		}
	}
}
