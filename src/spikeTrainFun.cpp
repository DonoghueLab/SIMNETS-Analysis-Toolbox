/**
 * @file spikeTrainFun.cpp
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
 */
#ifdef MATLAB_MEX_FILE
#include "mex.h"
#else
#endif
#include <vector>
#include <iostream>
#include <fstream>
#include <limits>
#include <array>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <iterator>
#include <functional>

/*
// Use this code to time execution:
#include <chrono>
std::chrono::high_resolution_clock::time_point t1, t2;
t1 = std::chrono::high_resolution_clock::now();
t2 = std::chrono::high_resolution_clock::now();
std::cout << "\tTime: "<< std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count() <<"us.\n";
*/

//using namespace std;
#include "spikeTrainFun.hpp"

template<typename T>
std::ostream& operator<< (std::ostream& out, const std::vector<T>& v) {
	out << "[";
	size_t last = v.size() - 1;
	for(size_t i = 0; i < v.size(); ++i) {
		out << v[i];
		if (i != last)
			out << ", ";
	}
	out << "]";
	return out;
}


SpikeTrain::SpikeTrain ()
{
	t_ = std::make_shared<SpkData_t>();
	win_t_[0] = - std::numeric_limits<SpkDatum_t>::infinity();
	win_t_[1] = std::numeric_limits<SpkDatum_t>::infinity();
	applyWin();
};

SpikeTrain::SpikeTrain (const SpkData_t& t)
{
	t_ = std::make_shared<SpkData_t>(t);
	win_t_[0] = - std::numeric_limits<SpkDatum_t>::infinity();
	win_t_[1] = std::numeric_limits<SpkDatum_t>::infinity();
	applyWin();
}

SpikeTrain::SpikeTrain (const SpkData_t& t, const std::array<SpkDatum_t, 2> win_t) : win_t_(win_t)
{
	t_ = std::make_shared<SpkData_t>(t);
	applyWin();
}

SpikeTrain::SpikeTrain (const SpkData_t& t, const std::array<SpkDatum_t, 2> win_t, const SpkDatum_t offset) : win_t_(win_t), _offset(offset)
{
	t_ = std::make_shared<SpkData_t>(t);
	applyWin();
}

SpikeTrain::SpikeTrain (const std::shared_ptr<SpkData_t> t_ptr) : t_(t_ptr)
{
	win_t_[0] = - std::numeric_limits<SpkDatum_t>::infinity();
	win_t_[1] = std::numeric_limits<SpkDatum_t>::infinity();
	applyWin();
}

SpikeTrain::SpikeTrain (const std::shared_ptr<SpkData_t> t_ptr, const std::array<SpkDatum_t, 2> win_t): t_(t_ptr), win_t_(win_t)
{
	applyWin();
}

SpikeTrain::SpikeTrain (const std::shared_ptr<SpkData_t> t_ptr, const std::array<SpkDatum_t, 2> win_t, const SpkDatum_t offset): t_(t_ptr), win_t_(win_t), _offset(offset)
{
	applyWin();
}

SpikeTrain::SpikeTrain (const std::shared_ptr<SpkData_t> t_ptr, const SpkDatum_t winStart, const SpkDatum_t winLength) :t_(t_ptr)
{
	std::array<SpkDatum_t, 2> win = {{winStart, winStart + winLength}};
	setWin(win);
}

SpikeTrain::SpikeTrain (const std::shared_ptr<SpkData_t> t_ptr, const SpkDatum_t winStart, const SpkDatum_t winLength, const SpkDatum_t offset) :t_(t_ptr), _offset(offset)
{
	std::array<SpkDatum_t, 2> win = {{winStart, winStart + winLength}};
	setWin(win);
}

SpikeTrain::SpikeTrain (const SpikeTrain& other)
{
	t_ = other.t_;
	win_t_ = other.win_t_;
	applyWin();
}

SpkData_t SpikeTrain::spikesInWindow() const
{
	SpkData_t w(win_i_[0], win_i_[1] );
	SpkDatum_t offset = _offset;
	if (! std::isinf(win_t_[0])) {
		offset -= win_t_[0];
	}
	for (auto &s : w) {
		s += offset;
	}
	return SpkData_t(w);
}

void SpikeTrain::print() const
{
	print(std::cout);
}

void SpikeTrain::print(std::ostream& os) const
{
	os << "Spk (" << length() <<"): " << std::setiosflags(std::ios::fixed) << std::setprecision(3);
	os << "Win [" << win_t_[0] << ", " << win_t_[1] << "]; ";
	SpkDatum_t subt = (std::isinf(win_t_[0])) ? 0 : win_t_[0];

	for (SpkData_t::iterator i = win_i_[0]; i != win_i_[1]; ++ i) {
		if (i != win_i_[0]) {
			os << " ";
		}
		os << (*i - subt);
	}
	os << "\n";
}

double SpikeTrain::distVPTo(const SpikeTrain& otherSt, const double q) const
{
	if (this == &otherSt) {
		return 0.0;
	}
	std::vector<double>::size_type n1 = length();
	std::vector<double>::size_type n2 = otherSt.length();
	std::vector<double>::size_type i, j, ind;
	double d;
	SpkDatum_t subt1 = (std::isinf(win_t_[0])) ? 0 : win_t_[0];
	SpkDatum_t subt2 = (std::isinf(otherSt.win_t_[0])) ? 0 : otherSt.win_t_[0];


	if (n1 == 0) { return (double)n2;}
	if (n2 == 0) { return (double)n1;}

	std::vector<double> dMat((n1 + 1) * (n2 + 1));
	for (i = 1; i <= n1; ++ i) {
		dMat[i] = (double) i;
	}
	for (i = 1; i <= n2; ++ i) {
		dMat[i * (n1 + 1)] = (double) i;
	}
	for (i = 1; i <= n1; ++i) {
		for (j = 1; j <= n2; ++j) {
			d = fabs(*(win_i_[0] + i - 1) - subt1 - *(otherSt.win_i_[0] + j -1) + subt2);
			ind = i + j * (n1 + 1);
			dMat[ind] = std::min(dMat[ind - (n1 + 1)] + 1.0, std::min(dMat[ind - 1] + 1.0, dMat[ind - 1 - (n1 + 1)] + q*d ));
		}
	}
    d = dMat[n1 * n2 + n1 + n2];
	return d;
}

void SpikeTrain::applyWin()
{
	// If no spikes, make window iterators point to nowhere
	if (t_->size() == 0) {
		win_i_[0] = t_->end();
		win_i_[1] = t_->end();
		return;
	}
	// Try to find first spike in window
	win_i_[0] = std::find_if(t_->begin(), t_->end(), [this](const double s) {return this->sInWin(s);});
	if (win_i_[0] != t_->end()) {
		// If found, find last one. We start at the point we first found, then
		// move to the first spike that falls out of the window. That should be
		// way faster than searching from end, if underlying spike trains are
		// much longer than typical windows. Leave previous code for comparison
		win_i_[1] = std::find_if(win_i_[0], t_->end(), [this](const double s) {return ! this->sInWin(s);});
		//auto a = std::find_if(t_->rbegin(), t_->rend(), [this](const double s) {return this->sInWin(s);});
		//win_i_[1] = (a).base();
	}
	else {
		// no spike in window, so iterators point to nowhere
		win_i_[0] = t_->end();
		win_i_[1] = t_->end();
		return;
	}
}

void SpikeTrain::setWin(const std::array<double, 2>& win)
{
	win_t_ = win;
	applyWin();
}

void SpikeTrain::setWin(const double winStart, const double winLength)
{
	std::array<double, 2> win = {{winStart, winStart + winLength}};
	setWin(win);
}

std::ostream& operator<<(std::ostream& os, const SpikeTrain& st)
{
	st.print(os);
	return os;
}

SpikeTrainList_t getSpikeTrainsFromDataWithWindow(const SpkDataList_t& tl, const std::array<SpkDatum_t, 2> win)
{
	SpikeTrainList_t::size_type nSt = tl.size();
	SpikeTrainList_t sptl;
	sptl.resize(nSt);
	for (SpikeTrainList_t::size_type i = 0; i < nSt; ++i) {
		sptl[i] = SpikeTrain(tl[i], win);
	}
	return sptl;
}

SpikeTrainList_t getSpikeTrainsFromDataWithWindow(const SpkDataList_t& tl, const SpkDatum_t winStart, const SpkDatum_t winLength)
{
	std::array<double, 2> win = {{winStart, winStart + winLength}};
	return getSpikeTrainsFromDataWithWindow(tl, win);
}

SpikeTrainList_t getSpikeTrainsFromDataWithWindow(const SpkData_t& t, const SpkData_t& aTimes, const std::array<SpkDatum_t, 2> win)
{
	return getSpikeTrainsFromDataWithWindow(std::make_shared<SpkData_t>(t), aTimes, win);
}

SpikeTrainList_t getSpikeTrainsFromDataWithWindow(const SpkData_t& t, const SpkData_t& aTimes, const SpkDatum_t winStart, const SpkDatum_t winLength)
{
	std::array<double, 2> win = {{winStart, winStart + winLength}};
	return getSpikeTrainsFromDataWithWindow(t, aTimes, win);
}

std::vector<SpikeTrainList_t> getSpikeTrainsFromDataWithWindow(const SpkDataList_t& tl, const SpkData_t& aTimes, const SpkDatum_t winStart, const SpkDatum_t winLength)
{
	std::vector<SpikeTrainList_t> outV(tl.size());
	for	(size_t i = 0; i < tl.size(); ++i) {
		outV[i] = getSpikeTrainsFromDataWithWindow(tl[i], aTimes, winStart, winLength);
	}
	return outV;
}


SpikeTrainList_t getSpikeTrainsFromDataWithWindow(const std::shared_ptr<SpkData_t>& t, const SpkData_t& aTimes, const std::array<SpkDatum_t, 2> win)
{
	SpikeTrainList_t::size_type nSt = aTimes.size();
	SpikeTrainList_t sptl;
	sptl.resize(nSt);

	#pragma omp parallel for shared(sptl)
	for (SpkData_t::size_type i = 0; i < nSt; ++i) {
		std::array<SpkDatum_t, 2> this_win;
		this_win[0] = win[0] + aTimes[i];
		this_win[1] = win[1] + aTimes[i];
		sptl[i] = SpikeTrain(t, this_win, win[0]);
	}
	return sptl;
}

SpikeTrainList_t getSpikeTrainsFromDataWithWindow(const std::shared_ptr<SpkData_t>& t, const SpkData_t& aTimes, const SpkDatum_t winStart, const SpkDatum_t winLength)
{
	std::array<double, 2> win = {{winStart, winStart + winLength}};
	return getSpikeTrainsFromDataWithWindow(t, aTimes, win);
}


SpikeTrainList_t getSpikeTrainsFromDataWithWindow(const SpkDataPtrList_t& tl, const std::array<SpkDatum_t, 2> win)
{
	SpkDataPtrList_t::size_type nSt = tl.size();
	SpikeTrainList_t sptl;
	sptl.resize(nSt);
	#pragma omp parallel for shared(sptl)
	for (SpkDataPtrList_t::size_type i = 0; i < nSt; ++i) {
		sptl[i] = SpikeTrain(tl[i], win);
	}
	return sptl;
}

SpikeTrainList_t getSpikeTrainsFromDataWithWindow(const SpkDataPtrList_t& tl, const SpkDatum_t winStart, const SpkDatum_t winLength)
{
	std::array<double, 2> win = {{winStart, winStart + winLength}};
	return getSpikeTrainsFromDataWithWindow(tl, win);
}


std::vector<SpikeTrainList_t> getSpikeTrainsFromDataWithWindow(const SpkDataPtrList_t& tl, const SpkData_t& aTimes, const SpkDatum_t winStart, const SpkDatum_t winLength)
{
	std::vector<SpikeTrainList_t> outV(tl.size());
	#pragma omp parallel for shared(outV)
	for	(size_t i = 0; i < tl.size(); ++i) {
		outV[i] = getSpikeTrainsFromDataWithWindow(tl[i], aTimes, winStart, winLength);
	}
	return outV;
}


std::vector<std::vector<std::vector<size_t>>> getSpikeDataEventWindowCounts(const SpkDataList_t &spkd, const SpkDataList_t &events, const SpkData_t &offsets, const SpkData_t &lengths)
{
	const size_t n_spkt = spkd.size(), n_events = events.size();
	std::vector<std::vector<std::vector<size_t>>> counts;
	counts.resize(n_spkt);

	#pragma omp parallel for shared(counts)
	for (size_t i = 0; i < n_spkt; ++i) {
		counts[i].resize(n_events);
		for (size_t j = 0; j < n_events; ++j) {
			// Generate spike trains for each spike train and event list combination, and windows
			SpikeTrainList_t spkTL = getSpikeTrainsFromDataWithWindow(spkd[i], events[j], offsets[j], lengths[j]);

			counts[i][j].resize(events[j].size());
			std::transform(spkTL.begin(), spkTL.end(), counts[i][j].begin(), std::bind(&SpikeTrain::size, std::placeholders::_1));
		}
	}
	return counts;
}

std::vector<std::vector<SpkDataList_t>> getSpikeDataInEventWindows(const SpkDataList_t &spkd, const SpkDataList_t &events, const SpkData_t &offsets, const SpkData_t &lengths)
{
	const size_t n_spkt = spkd.size(), n_events = events.size();
	std::vector<std::vector<SpkDataList_t>> spikeWindows;
	spikeWindows.resize(n_spkt);

	#pragma omp parallel for shared(spikeWindows)
	for (size_t i = 0; i < n_spkt; ++i) {
		spikeWindows[i].resize(n_events);
		for (size_t j = 0; j < n_events; ++j) {
			// Generate spike trains for each spike train and event list combination, and windows
			SpikeTrainList_t spkTL = getSpikeTrainsFromDataWithWindow(spkd[i], events[j], offsets[j], lengths[j]);

			spikeWindows[i][j].resize(events[j].size());
			std::transform(spkTL.begin(), spkTL.end(), spikeWindows[i][j].begin(), std::bind(&SpikeTrain::spikesInWindow, std::placeholders::_1));
		}
	}
	return spikeWindows;
}

Dist_t getDistMat(const SpikeTrainList_t& spikeTrains, const double q)
{
	Dist_t::size_type nSt = spikeTrains.size();
	Dist_t::size_type nDmat = nSt * (nSt - 1) / 2;
	Dist_t::size_type k = 0;
	Dist_t dmat(nDmat);
	//#pragma omp parallel for
	for (Dist_t::size_type i = 0; i < nSt - 1; ++ i) {
		for (Dist_t::size_type j = i + 1; j < nSt; ++ j) {
			dmat[k] = spikeTrains[i].distVPTo(spikeTrains[j], q);
			++ k;
		}
	}
	return dmat;
}

Dist_t::iterator getDistMat(const SpikeTrainList_t& spikeTrains, const double q,
	Dist_t::iterator dBegin, const Dist_t::const_iterator dEnd)
{
	Dist_t::size_type nSt = spikeTrains.size();
	//#pragma omp parallel for
	for (Dist_t::size_type i = 0; i < nSt - 1; ++ i) {
		for (Dist_t::size_type j = i + 1; j < nSt; ++ j) {
			if (dBegin != dEnd) {
				*dBegin = spikeTrains[i].distVPTo(spikeTrains[j], q);
				++ dBegin;
			}
		}
	}
	return dBegin;
}

arma::mat getFullDistMat(const SpikeTrainList_t& spikeTrains, const double q)
{
	size_t nSt = spikeTrains.size();
	arma::mat dmat(nSt, nSt);
	for (size_t i = 0; i < nSt; ++ i) {
		dmat(i, i) = 0.0;
	}
	#pragma omp parallel for shared(dmat)
	for (size_t i = 0; i < nSt - 1; ++ i) {
		for (size_t j = i + 1; j < nSt; ++ j) {
			dmat(i, j) = spikeTrains[i].distVPTo(spikeTrains[j], q);
			dmat(j, i) = dmat(i, j);
		}
	}
	return dmat;
}


Dist_t getDistMat(const SpikeTrainList_t& spikeTrains1, const SpikeTrainList_t& spikeTrains2, const double q)
{
	Dist_t::size_type nSt1 = spikeTrains1.size();
	Dist_t::size_type nSt2 = spikeTrains2.size();

	Dist_t::size_type nDmat = nSt1 * nSt2;
	Dist_t dmat(nDmat);
	#pragma omp parallel for shared(dmat)
	for (Dist_t::size_type j = 0; j < nSt2; ++ j) {
		for (Dist_t::size_type i = 0; i < nSt1; ++ i) {
			dmat[i + j * nSt1] = spikeTrains1[i].distVPTo(spikeTrains2[j], q);
		}
	}
	return dmat;
}

Dist_t::iterator getDistMat(const SpikeTrainList_t& spikeTrains1, const SpikeTrainList_t& spikeTrains2,
	const double q, Dist_t::iterator dBegin, const Dist_t::const_iterator dEnd)
{
	Dist_t::size_type nSt1 = spikeTrains1.size();
	Dist_t::size_type nSt2 = spikeTrains2.size();

	for (Dist_t::size_type j = 0; j < nSt2; ++ j) {
		for (Dist_t::size_type i = 0; i < nSt1; ++ i) {
			if (dBegin == dEnd) { return dBegin;}
			*dBegin = spikeTrains1[i].distVPTo(spikeTrains2[j], q);
			++ dBegin;
		}
	}
	return dBegin;
}

arma::mat getFullDistMat(const SpikeTrainList_t& spikeTrains1, const SpikeTrainList_t& spikeTrains2, const double q)
{
	size_t nSt1 = spikeTrains1.size();
	size_t nSt2 = spikeTrains2.size();
	arma::mat dmat(nSt1, nSt2);

	for (size_t j = 0; j < nSt2; ++ j) {
		for (size_t i = 0; i < nSt1; ++ i) {
			dmat(i, j) = spikeTrains1[i].distVPTo(spikeTrains2[j], q);
		}
	}
	return dmat;
}


// baseDmat should be of size: aTimes.size() x [tl.size() * (tl.size() - 1) / 2]
// baseSpikes should be of size aTimes.size() x tl.size()
void getSSIMDMatBetweenNeurons(const SpkDataList_t& tl, const SpkData_t& aTimes,
	const SpkDatum_t winStart, const SpkDatum_t winLength, const double q,
	Dist_t::iterator dBegin, const Dist_t::const_iterator dEnd,
	SpkDataList_t::iterator sBegin, const SpkDataList_t::const_iterator sEnd)
{
	Dist_t::iterator distPos = dBegin;
	SpkDataList_t::iterator spktPos = sBegin;

	for (auto &t: aTimes) {
		SpikeTrainList_t st = getSpikeTrainsFromDataWithWindow(tl, t + winStart, winLength);
		distPos = getDistMat(st, q, distPos, dEnd);
		for (auto &s: st) {
			if (spktPos != sEnd) {
				*spktPos = s.spikesInWindow();
				++spktPos;
			}
		}
	}
}

arma::mat getSSIMDMatBetweenNeurons(const SpkDataList_t& tl, const SpkData_t& aTimes,
	const SpkDatum_t winStart, const SpkDatum_t winLength, const double q,
	SpkDataList_t::iterator sBegin, const SpkDataList_t::const_iterator sEnd)
{
	const size_t n_tp = aTimes.size();
	const size_t n_u = tl.size();
	SpkDataList_t::iterator spktPos = sBegin;

	arma::mat dmat(n_u, n_tp * n_u);
	for (size_t i = 0; i < aTimes.size(); ++i) {
		SpikeTrainList_t st = getSpikeTrainsFromDataWithWindow(tl, aTimes[i] + winStart, winLength);
		dmat.submat(0, i * n_u, n_u - 1, (i + 1) * n_u - 1) = getFullDistMat(st, q);
		for (auto &s: st) {
			if (spktPos != sEnd) {
				*spktPos = s.spikesInWindow();
				++spktPos;
			}
		}
	}
	return dmat;
}

arma::mat getSSIMDMatBetweenNeurons(const SpkDataPtrList_t& tl, const SpkData_t& aTimes,
	const SpkDatum_t winStart, const SpkDatum_t winLength, const double q,
	SpkDataList_t::iterator sBegin, const SpkDataList_t::const_iterator sEnd)
{
	const size_t n_tp = aTimes.size();
	const size_t n_u = tl.size();
	SpkDataList_t::iterator spktPos = sBegin;

	arma::mat dmat(n_u, n_tp * n_u);
	for (size_t i = 0; i < aTimes.size(); ++i) {
		SpikeTrainList_t st = getSpikeTrainsFromDataWithWindow(tl, aTimes[i] + winStart, winLength);
		dmat.submat(0, i * n_u, n_u - 1, (i + 1) * n_u - 1) = getFullDistMat(st, q);
		for (auto &s: st) {
			if (spktPos != sEnd) {
				*spktPos = s.spikesInWindow();
				++spktPos;
			}
		}
	}
	return dmat;
}



// baseDmat should be of size: tl.size() x [aTimes.size() * (aTimes.size() - 1) / 2]
// baseSpikes should be of size tl.size() x aTimes.size()
void getSSIMDMatBetweenTimePoints(const SpkDataList_t& tl, const SpkData_t& timePoints,
	const SpkDatum_t winStart, const SpkDatum_t winLength, const double q,
	Dist_t::iterator dBegin, const Dist_t::const_iterator dEnd,
	SpkDataList_t::iterator sBegin, const SpkDataList_t::const_iterator sEnd)
{
	SpikeTrainList_t st;
	Dist_t::iterator distPos = dBegin;
	SpkDataList_t::iterator spktPos = sBegin;
	for (auto &u: tl) {
		st = getSpikeTrainsFromDataWithWindow(u, timePoints, winStart, winLength);
		distPos = getDistMat(st, q, distPos, dEnd);
		for (auto &s: st) {
			if (spktPos != sEnd) {
				*spktPos = s.spikesInWindow();
				++spktPos;
			}
		}
	}
}


arma::mat getSSIMDMatBetweenTimePoints(const SpkDataList_t& tl, const SpkData_t& timePoints,
const SpkDatum_t winStart, const SpkDatum_t winLength, const double q,
SpkDataList_t::iterator sBegin, const SpkDataList_t::const_iterator sEnd)
{
	SpkDataList_t::iterator spktPos = sBegin;
	size_t i = 0;
	const size_t n_tp = timePoints.size();
	const size_t n_u = tl.size();
	arma::mat dmat(n_tp, n_tp * n_u);

	for (auto &u: tl) {
		SpikeTrainList_t st = getSpikeTrainsFromDataWithWindow(u, timePoints, winStart, winLength);
		dmat.submat(0, i * n_tp, n_tp - 1, (i + 1) * n_tp - 1) = getFullDistMat(st, q);
		for (auto &s: st) {
			if (spktPos != sEnd) {
				*spktPos = s.spikesInWindow();
				++spktPos;
			}
		}
		++i;
	}
	return dmat;
}

arma::mat getSSIMDMatBetweenTimePoints(const SpkDataPtrList_t& tl, const SpkData_t& timePoints,
const SpkDatum_t winStart, const SpkDatum_t winLength, const double q,
SpkDataList_t::iterator sBegin, const SpkDataList_t::const_iterator sEnd)
{
	SpkDataList_t::iterator spktPos = sBegin;
	size_t i = 0;
	const size_t n_tp = timePoints.size();
	const size_t n_u = tl.size();
	arma::mat dmat(n_tp, n_tp * n_u);

	for (auto &u: tl) {
		SpikeTrainList_t st = getSpikeTrainsFromDataWithWindow(u, timePoints, winStart, winLength);
		dmat.submat(0, i * n_tp, n_tp - 1, (i + 1) * n_tp - 1) = getFullDistMat(st, q);
		for (auto &s: st) {
			if (spktPos != sEnd) {
				*spktPos = s.spikesInWindow();
				++spktPos;
			}
		}
		++i;
	}
	return dmat;
}

// function to copy lower triag matrix to square matrix
// dest must be allocated and large enough to hold all distance data
size_t copyPDistToSquare(const Dist_t &dm, double *dest)
{
	// determine number of rows & columns
	size_t nIt = (1 + sqrt(1 + 8*dm.size()))/2;
	size_t ind = 0;
	if (dest != NULL) {
		// set diagonal to 0.0
		for (size_t i = 0; i < nIt; ++i) {
			dest[i + nIt * i] = 0.0;
		}
		ind = 0;
		// copy entries from dm to both upper and lower triangle in dest
		for (size_t i = 0; i < nIt -1; ++i) {
			for (size_t j = i + 1; j < nIt; ++j) {
				dest[j + nIt * i] = dm[ind];
				dest[i + nIt * j] = dm[ind];
				++ ind;
			}
		}
	}
	return ind;
}

arma::mat copyNTimesPDistToMat(const Dist_t &dm, const size_t n_reps, const size_t n_pw)
{
	arma::mat dmat;
	size_t ind, n_pwds = n_pw * (n_pw - 1) / 2;
	if (dm.size() != n_reps * n_pwds) {
		std::cout << "oops.\n";
		return dmat;
	}
	dmat.set_size(n_pw, n_reps * n_pw );
	for (size_t r = 0; r < n_reps; ++r) {
		for (size_t i = 0; i < n_pw; ++i) {
			dmat(i, i + r*n_pw) = 0.0;
		}
		ind = 0;
		for (size_t i = 0; i < n_pw; ++i) {
			for (size_t j = i + 1; j < n_pw; ++j) {
				dmat(i, j + r*n_pw) = dm[ind + n_pwds * r];
				dmat(j, i + r*n_pw) = dm[ind + n_pwds * r];
				++ ind;
			}
		}
	}
	return dmat;
}

void sanitizeSpikeData(SpkDataList_t& tl)
{
	for (auto &t: tl) {
		sanitizeSpikeData(t);
	}
}

void sanitizeSpikeData(SpkData_t& t)
{
	std::sort(t.begin(), t.end());
}

void sanitizeSpikeData(SpkDataPtrList_t& tl)
{
	for (auto &t: tl) {
		sanitizeSpikeData(*t);
	}
}


void sanitizeSpikeData(std::shared_ptr<SpkData_t> t)
{
	sanitizeSpikeData(*t);
}



SpkDataPtrList_t loadSpikeTrains(const char* fn)
{
	uint64_t n, m;
	double t;
	size_t i, j;
	SpkDataPtrList_t spkl;

	SpkData_t spkt;
	std::ifstream infile;
	infile.open(fn, std::ios::binary | std::ios::in);
	if (infile.is_open()) {
		infile.read((char*)&n, 8);
		spkl.resize(n);
		std::cout << "We have "<<n<<" spike trains\n";
		for (i = 0; i < n; ++i) {
			infile.read((char*)&m, 8);
			spkt.resize(m);
			std::cout << "Spkt "<<i<<" has "<<m<<" spikes\n";
			for (j = 0; j < m; ++j) {
				infile.read((char*)&t, 8);
				spkt[j] = t;
			}
			spkl[i] = std::make_shared<SpkData_t>(spkt);
		}
	}
	infile.close();
	return spkl;
}
