/**
 * @file getSSIMSprojection.cpp
 * @brief MATLAB wrapper for SpikeTrain class, to apply tSNE projection on new data, (mostly) replicating getSSIMSprojection.m
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
#include <armadillo>
#include "armaMex.hpp"
#include "spikeTrainFun.hpp"

#define MAZ_SZ (~(mwSize)0)

/**	MATLAB wrapper for SpikeTrain, replicating getSSIMSprojection.m.
 * function [SSIMS_coordinates] = getSSIMSprojection(spikeTrains, evTimes, startTimes, winLen, q, tSNE_transform, basespiketrains)
 *
 * 	When compiled will produce a mex function which takes three arguments:
 * 	@param spikeTrains cell array of n vectors of doubles, representing spike trains from n neurons
 * 	@param evTimes vector of m doubles, representing times at which distances are taken
 * 	@param startTimes beginning of window relative to evTimes (double vector)
 * 	@param winLen window length (double scalar)
 * 	@param q Cost for moving a spike (temporal resolution; double scalar)
 * 	@param tSNE_transform tSNE_transform as returned by getSSIMS
 * 	@param basespiketrains Spike trains to be compared agains
 *
 * 	@return 	SSIMS_coordinates, a matrix of size length(evTimes) * size(tSNE_transform, 2) * length(startTimes)
**/
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	const mxArray *spikeTrains, *evTimes, *startTimes, *baseSpikeTrains, *tSNE_transform;  // input spike train matrices

	mwSize nTmpList;
	double q, winLen = 1.0, *startTimesDbl;
	bool isNotDouble, isNotCell;
	arma::mat aBaseData;
	mwIndex badIndex;
	// check for proper number of arguments
	if (nrhs < 7) {
		mexErrMsgIdAndTxt("SSIMSToolbox:getSSIMSprojection:nrhs", "This function expects at least 7 arguments.");
	}

	// Get input variables and make sure they are of correct type
	spikeTrains = prhs[0];
	if( !mxIsCell(spikeTrains) || mxIsComplex(spikeTrains)) {
		mexErrMsgIdAndTxt("SSIMSToolbox:getSSIMSprojection:notDouble", "Expects cell array of spike trains as argument 1.");
	}
	const mwSize nSpikeTrains = mxGetNumberOfElements(spikeTrains);
	isNotDouble = false;
	for (mwIndex i = 0; i < nSpikeTrains; ++i) {
		if (!mxIsDouble(mxGetCell(spikeTrains, i))) {
			isNotDouble = true;
			badIndex = i;
			break;
		}
	}
	if (isNotDouble) {
		mexErrMsgIdAndTxt("SSIMSToolbox:getSSIMSprojection:notDouble","Element %i of spike train cell has class %s. We expect double!", badIndex, mxGetClassName(mxGetCell(spikeTrains, badIndex)));
	}

	// make sure other arguments have correct format
	evTimes = prhs[1];
	if ( !mxIsDouble(evTimes) || mxIsComplex(evTimes)) {
		mexErrMsgIdAndTxt("SSIMSToolbox:getSSIMSprojection:notDouble2", "Times (argument 2) must be real doubles.");
	}

	startTimes = prhs[2];
	if ( !mxIsDouble(startTimes) || mxIsComplex(startTimes)) {
		mexErrMsgIdAndTxt("SSIMSToolbox:getSSIMSprojection:notDouble3","Window start times must be double.");
	}

	if ((nrhs > 3) && ( !mxIsDouble(prhs[3]) || mxIsComplex(prhs[3]) || mxGetNumberOfElements(prhs[3])!=1)) {
		mexErrMsgIdAndTxt("SSIMSToolbox:getSSIMSprojection:notScalar4","Window length must be scalar.");
	}
	winLen = mxGetScalar(prhs[3]);

	if ((nrhs > 4) && ( !mxIsDouble(prhs[4]) || mxIsComplex(prhs[4]) || mxGetNumberOfElements(prhs[4])!=1)) {
		mexErrMsgIdAndTxt("SSIMSToolbox:getSSIMSprojection:notScalar5","Cost q must be scalar.");
	}
	q = mxGetScalar(prhs[4]);

	tSNE_transform = prhs[5];
	if ((nrhs > 5) && ( !mxIsDouble(tSNE_transform) || mxIsComplex(tSNE_transform))) {
		mexErrMsgIdAndTxt("SSIMSToolbox:getSSIMSprojection:notDouble6","tSNE_transform must be double.");
	}

	baseSpikeTrains = prhs[6];
	if (!mxIsCell(baseSpikeTrains) || mxIsComplex(baseSpikeTrains)) {
		mexErrMsgIdAndTxt("SSIMSToolbox:getSSIMSprojection:notDouble", "Expects cell array of spike trains as argument 7.");
	}
	const mwSize nBaseSpikeTrains = mxGetNumberOfElements(baseSpikeTrains);
	if (nBaseSpikeTrains != nSpikeTrains) {
		mexErrMsgIdAndTxt("SSIMSToolbox:getSSIMSprojection:notDouble", "Spike train cell array and baseSpikeTrains have to have same number of elements (units).");
	}
	isNotDouble = false;
	isNotCell = false;
	for (mwIndex i = 0; i < nBaseSpikeTrains; ++i) {
		if (!mxIsCell(mxGetCell(baseSpikeTrains, i))) {
			badIndex = i;
			isNotCell = true;
			break;
		}
		else {
			for (mwIndex j = 0; j < mxGetNumberOfElements(mxGetCell(baseSpikeTrains, i)); ++j) {
				if (!mxIsDouble(mxGetCell(mxGetCell(baseSpikeTrains, i), j))) {
					badIndex = i;
					isNotDouble = true;
					break;
				}
			}
		}
	}
	if (isNotDouble) {
		mexErrMsgIdAndTxt("SSIMSToolbox:getSSIMSprojection:notDouble","An element of base spike train cell  %i has class %s. We expect double!", badIndex, mxGetClassName(mxGetCell(baseSpikeTrains, badIndex)));
	}
	if (isNotCell) {
		mexErrMsgIdAndTxt("SSIMSToolbox:getSSIMSprojection:notCell","Element %i of base spike train cell has class %s. We expect a cell!", badIndex, mxGetClassName(mxGetCell(baseSpikeTrains, badIndex)));
	}

	const mwSize nTimePoints = mxGetNumberOfElements(evTimes);
	const mwSize nStartTimes = mxGetNumberOfElements(startTimes);
	const mwSize nDims = mxGetN(tSNE_transform);
	SpkDataPtrList_t tl(nSpikeTrains);

	mxArray * tmpStl, *tmpSt;
	SpkData_t tmpSpkd;
	// Compile raw spike data
	for (SpkDataPtrList_t::iterator s = tl.begin(); s != tl.end(); ++s) {
		tmpSt = mxGetCell(spikeTrains, s - tl.begin());
		*s = std::make_shared<SpkData_t> (mxGetPr(tmpSt), mxGetPr(tmpSt) + mxGetNumberOfElements(tmpSt));
	}
	sanitizeSpikeData(tl);
	// Fill base spike trains
	std::vector<SpikeTrainList_t> baseSpikeTrainList(nBaseSpikeTrains);
	mwSize firstSize = MAZ_SZ;
	for (mwIndex i = 0; i < nBaseSpikeTrains; ++i) {
		tmpStl = mxGetCell(baseSpikeTrains, i);
		nTmpList = mxGetNumberOfElements(tmpStl);
		if (i == 0) {
			firstSize = nTmpList;
		}
		else if (firstSize != nTmpList) {
			mexErrMsgIdAndTxt("SSIMSToolbox:getSSIMSprojection:baseSpkSizeMismatch","Mismatch between number of elements in basespiketrains: the first element has %i trains, the %ith one has %i. They have to have the same number of elements.", firstSize, i + 1, nTmpList);
		}
		SpikeTrainList_t tmpList(nTmpList);
		for (mwIndex j = 0; j < nTmpList; ++j) {
			tmpSt = mxGetCell(tmpStl, j);
			tmpSpkd = SpkData_t (mxGetPr(tmpSt), mxGetPr(tmpSt) + mxGetNumberOfElements(tmpSt));
			tmpList[j] = SpikeTrain(tmpSpkd);
		}
		baseSpikeTrainList[i] = tmpList;
	}
	const mwSize nBaseEvents = nTmpList;

	if (nBaseSpikeTrains * nBaseEvents != mxGetM(tSNE_transform)) {
		mexErrMsgIdAndTxt("SSIMSToolbox:getSSIMSprojection:sizeMismatch","Size mismatch between number of units (%i), number of events (%i), and transformation matrix (%ix%i, expected %ix%i).", nBaseSpikeTrains, nBaseEvents, mxGetM(tSNE_transform), nDims, nTimePoints*nSpikeTrains, nDims);
	}


	startTimesDbl = mxGetPr(startTimes);

	arma::cube SSIMSproj(nTimePoints, nDims, nStartTimes, arma::fill::zeros);
	const arma::mat tSNE_transform_mat = armaGetPr(tSNE_transform);

	// For each start time, get SSIMS projection
	#pragma omp parallel for schedule(dynamic) shared(SSIMSproj)
	for (mwIndex i = 0; i < nStartTimes; ++i) {
		std::vector<SpikeTrainList_t> projSpkt = getSpikeTrainsFromDataWithWindow(tl, SpkData_t(mxGetPr(evTimes), mxGetPr(evTimes) + mxGetNumberOfElements(evTimes)), startTimesDbl[i], winLen);
		arma::mat tmpDmat(nTimePoints, nBaseSpikeTrains * nBaseEvents);
		for (mwIndex j = 0; j < nBaseSpikeTrains; ++j) {
			tmpDmat.submat(0, j * nBaseEvents, nTimePoints - 1, (j + 1) * nBaseEvents - 1) = getFullDistMat(projSpkt[j], baseSpikeTrainList[j], q);
		}
		SSIMSproj.slice(i) = tmpDmat * tSNE_transform_mat;
	}

	if (nlhs > 0) {
		plhs[0] = armaCreateMxMatrix(SSIMSproj.n_rows, SSIMSproj.n_cols, SSIMSproj.n_slices);
		armaSetCubeData(plhs[0], SSIMSproj);

	}
}
