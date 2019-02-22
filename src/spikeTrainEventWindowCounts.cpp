/**
 * @file spikeTrainEventWindowCounts.cpp
 * \brief MATLAB wrapper for getSpikeDataEventWindowCounts()
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


mxArray * getMexArray(const std::vector<double>& v){
	mxArray * mx = mxCreateDoubleMatrix(v.size(), 1, mxREAL);
	std::copy(v.begin(), v.end(), mxGetPr(mx));
	return mx;
}

mxArray * getMexArray(const std::vector<double>& v, SpkDatum_t l){
	mxArray * mx = mxCreateDoubleMatrix(v.size(), 1, mxREAL);
	std::copy(v.begin(), v.end(), mxGetPr(mx));
	for (size_t i = 0; i < v.size(); ++i) {
		*(mxGetPr(mx) + i) = *(mxGetPr(mx) + i) / l;
	}
	return mx;
}


mxArray * getMexArray(const std::vector<size_t>& v) {
	std::vector<double> v_double(v.begin(), v.end());
	return getMexArray(v_double);
}

mxArray * getMexArray(const std::vector<size_t>& v, SpkDatum_t l) {
	std::vector<double> v_double(v.begin(), v.end());
	return getMexArray(v_double, l);
}

/**
 * MATLAB wrapper for getSpikeDataEventWindowCounts().
 *
 * When compiled will produce a mex function which takes five arguments:
 * @param spkT is a single spike train (vector of spike times [double]), or a cell array of spike trains
 * @param events is a vector of event times (double) or a cell array of event time vectors
 * @param offsets is a single offset time (relative to event times) for spike train window, or it is a vector with numel(events) elements (only if events is cell)
 * @param lengths is a single window length (relative to offsets) for spike train window, or it is a vector with numel(events) elements (only if events is cell)
 * @param normalize: {0, [1]} If true (!=0), spike counts will be normalized by window lengths
 * @return cell array; numel(spkT) x numel(events) (assuming spkT and events being cell arrays)
**/
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	// check for proper number of arguments
	if (nrhs < 2) {
		mexErrMsgIdAndTxt("SSIMSToolbox:spikeTrainEventWindowCounts:nrhs", "This function expects at least 2 inputs.");
	}
	mwSize n_spkt, n_events;
	SpkDataList_t spkd, events;
	bool normalize_counts = true;

	if (mxIsCell(prhs[0])) {
		n_spkt = mxGetNumberOfElements(prhs[0]);
		spkd.resize(n_spkt);
		for (SpkDataList_t::iterator s = spkd.begin(); s != spkd.end(); ++s) {
			mxArray* tmpSt = mxGetCell(prhs[0], s - spkd.begin());
			if (!mxIsDouble(tmpSt) || mxIsComplex(tmpSt)) {
				mexErrMsgIdAndTxt("SSIMSToolbox:spikeTrainEventWindowCounts:notDouble", "Spike trains must be real doubles!");
			}
			*s = SpkData_t (mxGetPr(tmpSt), mxGetPr(tmpSt) + mxGetNumberOfElements(tmpSt));
		}
	}
	else {
		if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0])) {
			mexErrMsgIdAndTxt("SSIMSToolbox:spikeTrainEventWindowCounts:notDouble", "Spike times must be real double!");
		}
		n_spkt = 1;
		spkd.push_back(SpkData_t(mxGetPr(prhs[0]), mxGetPr(prhs[0]) + mxGetNumberOfElements(prhs[0])));
	}
	sanitizeSpikeData(spkd);

	if (mxIsCell(prhs[1])) {
		n_events = mxGetNumberOfElements(prhs[1]);
		events.resize(n_events);
		for (SpkDataList_t::iterator s = events.begin(); s != events.end(); ++s) {
			mxArray* tmpSt = mxGetCell(prhs[1], s - events.begin());
			if (!mxIsDouble(tmpSt) || mxIsComplex(tmpSt)) {
				mexErrMsgIdAndTxt("SSIMSToolbox:spikeTrainEventWindowCounts:notDouble", "Event time points must be double!");
			}
			*s = SpkData_t (mxGetPr(tmpSt), mxGetPr(tmpSt) + mxGetNumberOfElements(tmpSt));
		}
	}
	else {
		if (!mxIsDouble(prhs[1]) || mxIsComplex(prhs[1])) {
			mexErrMsgIdAndTxt("SSIMSToolbox:spikeTrainEventWindowCounts:notDouble", "Event time points must be double!");
		}
		n_events = 1;
		events.push_back(SpkData_t(mxGetPr(prhs[1]), mxGetPr(prhs[1]) + mxGetNumberOfElements(prhs[1])));
	}

	SpkData_t offsets(n_events, 0.0), lengths(n_events, 1.0);
	if (nrhs > 2) {
		if (mxGetNumberOfElements(prhs[2]) == 1) {
			for (SpkData_t::iterator s = offsets.begin(); s != offsets.end(); ++s) {
				*s = mxGetScalar(prhs[2]);
			}
		}
		else if (mxGetNumberOfElements(prhs[2]) == n_events) {
			for (SpkData_t::iterator s = offsets.begin(); s != offsets.end(); ++s) {
				*s = *(mxGetPr(prhs[2]) +  (s - offsets.begin())) ;
			}
		}
		else {
			mexErrMsgIdAndTxt("SSIMSToolbox:spikeTrainEventWindowCounts:offsetsCount", "Offsets must have same number of elements as events, or be a scalar!");
		}
	}
	if (nrhs > 3) {
		if (mxGetNumberOfElements(prhs[3]) == 1) {
			for (SpkData_t::iterator s = lengths.begin(); s != lengths.end(); ++s) {
				*s = mxGetScalar(prhs[3]);
			}
		}
		else if (mxGetNumberOfElements(prhs[3]) == n_events) {
			for (SpkData_t::iterator s = lengths.begin(); s != lengths.end(); ++s) {
				*s = *(mxGetPr(prhs[3]) +  (s - lengths.begin())) ;
			}
		}
		else {
			mexErrMsgIdAndTxt("SSIMSToolbox:spikeTrainEventWindowCounts:lengthsCount", "Lengths must have same number of elements as events, or be a scalar!");
		}
	}
	if (nrhs > 4) {
		normalize_counts = mxGetScalar(prhs[4]) != 0;
	}



	mxArray *counts = NULL;
	if(nlhs > 0) {
		counts = mxCreateCellMatrix(n_spkt, n_events);
		plhs[0] = counts;

		std::vector<std::vector<std::vector<size_t>>> v_counts;
		v_counts = getSpikeDataEventWindowCounts(spkd, events, offsets, lengths);
		for (size_t i = 0; i < n_spkt; ++i) {
			for (size_t j = 0; j < n_events; ++j) {
				mwIndex subs[2] = {i, j};
				mxArray *cellItem;
				if (normalize_counts) {
					cellItem = getMexArray(v_counts[i][j], lengths[j]);
				}
				else {
					cellItem = getMexArray(v_counts[i][j]);
				}

				mxSetCell(counts, mxCalcSingleSubscript(counts, 2, subs), cellItem);
			}
		}
	}

}
