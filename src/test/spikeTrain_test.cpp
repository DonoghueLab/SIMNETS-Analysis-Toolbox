/**
 * @file SSIMS_test.cpp
 * @author  Jonas B. Zimmermann <jonas_zimmermann@brown.edu>
 * @version 0.1
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
#undef NDEBUG
#include <cassert>
#include <iostream>
#include <iomanip>
#include <iterator>
#include <algorithm>
#include <chrono>
using namespace std::chrono;
#include <armadillo>

extern "C" {
	#include "knnHelpers.h"
}
#include "tSNERunnerV.hpp"
#include "spikeTrainFun.hpp"
#include "mat.h"
#include "example_data_helper.hpp"

template<typename T>
std::ostream& operator<< (std::ostream& out, const std::vector<T>& v) {
	out << "[";
	size_t last = v.size() - 1;
	for(size_t i = 0; i < v.size(); ++i) {
		out << v[i];
		if (i != last) {
			out << ", ";
		}
	}
	out << "]";
	return out;
}

int test_spike_windows(const SpkDataList_t& spks)
{
	size_t n = std::min(10lu, spks[0].size());
	std::cout << "First " << n << " spike times of train 1:" << std::endl << std::fixed << std::setprecision(3);
	for (size_t i = 0; i < n; ++i) {
		std::cout << std::setw(12) << spks[0][i];
	}
	std::cout << std::endl;

	SpikeTrain s1(spks[0]);
	std::cout << "Spike train 1 contains " << std::setprecision(0) << s1.length() << " spikes." << std::endl;
	assert(11077 == s1.length());

	std::array<SpkDatum_t, 2> win = {{0.1, 0.4}};
	SpikeTrain s2(spks[0], win, 0);
	SpikeTrain s3(spks[0], win, 0.1);

	std::cout << "Spike train 2 with window [0.1s, 0.4s): " << s2 << "  new spike data from window: " << s2.spikesInWindow() << std::endl;
	std::cout << "Spike train 3 with window [0.1s, 0.4s) and 0.1s offset: " << s3 << "  new spike data from window: " << s3.spikesInWindow() << std::endl;
	std::cout << std::endl;


	s1.setWin(0, 3);

	std::cout << "Spike train window set to [0, 3)s. " << s1;

	SpkData_t event_times = {{0, 1.0, 1.5, 2.0, 2.1}};
	SpikeTrainList_t stlist;
	stlist = getSpikeTrainsFromDataWithWindow(s1.spikesInWindow(), event_times, -0.1, 0.5);
	std::cout << "We now apply times in window [-.1, .4) around events " << event_times << "." << std::endl << stlist;

	return 0;
}

int test_spike_windows_events(const SpkDataList_t& spks, const SpkDataList_t& cue_list)
{
	std::vector<SpikeTrainList_t> spkTrainsMov;
	spkTrainsMov = getSpikeTrainsFromDataWithWindow(spks, cue_list[2], -.2, 1);

	std::cout << "Spike trains extracted: " << spkTrainsMov.size() << std::endl;


	return 0;

}



int test_data(const char *file)
{
	std::cout << "Reading data from file " << file << std::endl;
	SpkDataList_t cue_list(3);
	SpkData_t trial_group;
	SpkDataList_t spikes;

	load_data(file, cue_list, trial_group, spikes);

	print_sample_structure(cue_list, trial_group, spikes);
	test_spike_windows(spikes);
	test_spike_windows_events(spikes, cue_list);


	return 0;
}

int main (int argc, char const *argv[])
{
	if (argc > 1) {
		return test_data(argv[1]);
	}
	else {
		return test_data("../examples/SSIMS_demo_data_center_out.mat");
	}
	return 0;
}
