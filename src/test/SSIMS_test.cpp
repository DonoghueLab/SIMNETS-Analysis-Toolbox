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

bool dsort (double x, double y) { return (x<y); }

int test_spike_tSNE_move(const SpkDataList_t& spks, const SpkDataList_t& cue_list, const SpkData_t& trial_group)
{
	const size_t nBaseSpikeTrains = spks.size();
	const size_t nEvents = trial_group.size();
	SpkDataList_t baseSpikesCo(nBaseSpikeTrains * nEvents);
	arma::mat dmatCo;
	double myQ = 10.0, kld, perplexity = 30.0, winLen = 1.0, winOffset = -0.1;
	size_t nDims = 3, nDimsClass = 10;

	high_resolution_clock::time_point t1, t2;

	t1 = high_resolution_clock::now();
	dmatCo = getSSIMDMatBetweenTimePoints(spks, cue_list[2], winOffset, winLen, myQ, baseSpikesCo.begin(), baseSpikesCo.end());
	t2 = high_resolution_clock::now();
	std::cout << "Calculating base distance matrix took "<< std::chrono::duration_cast<std::chrono::milliseconds>( t2 - t1 ).count() <<" ms.\n";

	std::cout<< "dmat size: "<< dmatCo.n_rows << " × "<< dmatCo.n_cols << "\n";

	t1 = high_resolution_clock::now();
	TSNERunnerV T1(dmatCo, nDims, perplexity);
	kld = T1.run();
	t2 = high_resolution_clock::now();
	std::cout << "Transformed data:\n" << T1.getTSNEData().rows(0, 3) << "...\n";
	std::cout << "Perplexity: " << std::fixed << std::setprecision(1) << T1.getPerplexity() << "\nKL d: " << std::setprecision(3) << kld << std::endl;
	std::cout << "3-dim tSNE took: " << std::setprecision(0) << std::chrono::duration_cast<std::chrono::milliseconds>( t2 - t1 ).count() << "ms\n";

	std::cout << "Making base spike trains ... ";
	std::vector<SpikeTrainList_t> basespiketrains(nBaseSpikeTrains);
	for (size_t i = 0; i < nBaseSpikeTrains; ++i) {
		basespiketrains[i].reserve(nEvents);
		for (size_t j = 0; j < nEvents; ++j) {
			basespiketrains[i].push_back(SpikeTrain(baseSpikesCo[i * nEvents + j]));
		}
	}
	std::cout << "Done." << std::endl;

	size_t nStartTimes = 41;
	arma::vec startTimes = arma::linspace(-1, 1, nStartTimes); // Make projection times -1:0.05:1

	arma::cube SSIMSproj(nEvents, nDims, nStartTimes, arma::fill::zeros);

	// For each start time, get SSIMS projection
 	std::cout << "Calculate projections: "<< std::endl;
	t1 = high_resolution_clock::now();
	#pragma omp parallel for schedule(dynamic) shared(SSIMSproj)
	for (size_t i = 0; i < nStartTimes; ++i) {
		std::vector<SpikeTrainList_t> projSpkt = getSpikeTrainsFromDataWithWindow(spks, cue_list[2], startTimes[i], winLen);
		arma::mat tmpDmat(nEvents, nBaseSpikeTrains * nEvents);
		for (size_t j = 0; j < nBaseSpikeTrains; ++j) {
			tmpDmat.submat(0, j * nEvents, nEvents - 1, (j + 1) * nEvents - 1) = getFullDistMat(projSpkt[j], basespiketrains[j], myQ);
		}
		SSIMSproj.slice(i) = T1.projectDataIntoTSNESpace(tmpDmat);
	}
	t2 = high_resolution_clock::now();
	std::cout << "Done. Projection of " << nStartTimes << " time steps took: " << std::chrono::duration_cast<std::chrono::milliseconds>( t2 - t1 ).count() << "ms\n";

	size_t i = 0;
	std::cout << "Projection at time " << std::setprecision(3) << startTimes[i] << ":\n" << SSIMSproj.slice(i).rows(0, 3)<< " ... " << std::endl << std::endl ;



	t1 = high_resolution_clock::now();
	TSNERunnerV T2(dmatCo, nDimsClass, perplexity);
	kld = T2.run();
	t2 = high_resolution_clock::now();
	std::cout << "Perplexity: " << std::setprecision(1) << T2.getPerplexity() << "\nKL d: " << std::setprecision(3) << kld << std::endl;
	std::cout << "10-dim tSNE took: " << std::setprecision(0) << std::chrono::duration_cast<std::chrono::milliseconds>( t2 - t1 ).count() << "ms\n";


	// Some (ugly, because C style) code to get K-NN classification
	const size_t k = 1;
	const size_t n_perm = 1000;

	std::vector<double> pperm(n_perm);
	std::vector<long long> classes(nEvents);
	for (size_t i = 0; i < nEvents; ++i) {
		classes[i] = trial_group[i];
	}
	std::vector<long long> knnclass(nEvents);
	std::vector<long long> uclasses(nEvents);
	ssize_t n_unique_classes = findUniqueItemsLL(&classes[0], nEvents, &uclasses[0]);
	uclasses.resize(n_unique_classes);

	std::vector<double> vm(nEvents * n_unique_classes);

	double perc;

	std::vector<double> dist_mat_tril(nEvents * (nEvents - 1) / 2);

	const arma::mat tSNEData = T2.getTSNEData();
	const double* tSNEData_mem = tSNEData.memptr();
	size_t n_cols = tSNEData.n_cols;

	Rpdist(tSNEData_mem, nEvents, n_cols, &dist_mat_tril[0]);

	srand((unsigned) time(NULL));
	perc = doKnnClassification(&dist_mat_tril[0], nEvents, &classes[0], k, n_perm, &knnclass[0], &vm[0], &pperm[0]);

	std::sort (pperm.begin(), pperm.end(), dsort);

	size_t eixl = (n_perm - 1) * 0.025;
	size_t eixm = (n_perm - 1) * 0.5;
	size_t eixh = (n_perm - 1) * 0.975;

	std::cout << "KNN classification percentage: " << std::setprecision(1) << perc << "%, with error bounds: [" << pperm[eixl] << ", " << pperm[eixm] << ", " << pperm[eixh] << "]." << std::endl;

	return 0;
}

int test_spike_tSNE_combined(const SpkDataList_t& spks, const SpkDataList_t& cue_list, const SpkData_t& trial_group)
{
	const size_t nBaseSpikeTrains = spks.size();
	const size_t nTrials = trial_group.size();

	arma::mat dmatCo;
	double myQ = 10.0, kld, perplexity = 30.0, winLen = 1.0, winOffset = -0.1;
	size_t nDims = 5;

	high_resolution_clock::time_point t1, t2;

	SpkData_t ev_t(cue_list[2]);
	for (SpkData_t::iterator it = ev_t.begin() ; it != ev_t.end(); ++it) {
		*it = *it + winOffset;
	}
	ev_t.insert(ev_t.end(), cue_list[0].begin(), cue_list[0].end());
	size_t nEvents = ev_t.size();
	SpkDataList_t baseSpikesCo(nBaseSpikeTrains * nEvents);

	t1 = high_resolution_clock::now();
	dmatCo = getSSIMDMatBetweenTimePoints(spks, ev_t, 0.0, winLen, myQ, baseSpikesCo.begin(), baseSpikesCo.end());
	t2 = high_resolution_clock::now();
	std::cout << "Calculating base distance matrix took "<< std::chrono::duration_cast<std::chrono::milliseconds>( t2 - t1 ).count() <<" ms.\n";

	std::cout<< "dmat size: "<< dmatCo.n_rows << " × "<< dmatCo.n_cols << "\n";

	t1 = high_resolution_clock::now();
	TSNERunnerV T1(dmatCo, nDims, perplexity);
	kld = T1.run();
	t2 = high_resolution_clock::now();
	std::cout << "Transformed data:\n" << T1.getTSNEData().rows(0, 3) << "...\n";
	std::cout << "Perplexity: " << std::fixed << std::setprecision(1) << T1.getPerplexity() << "\nKL d: " << std::setprecision(3) << kld << std::endl;
	std::cout << nDims << "-dim tSNE took: " << std::setprecision(0) << std::chrono::duration_cast<std::chrono::milliseconds>( t2 - t1 ).count() << "ms\n";


	std::cout << "Making base spike trains ... ";
	std::vector<SpikeTrainList_t> basespiketrains(nBaseSpikeTrains);
	for (size_t i = 0; i < nBaseSpikeTrains; ++i) {
		basespiketrains[i].reserve(nEvents);
		for (size_t j = 0; j < nEvents; ++j) {
			basespiketrains[i].push_back(SpikeTrain(baseSpikesCo[i * nEvents + j]));
		}
	}
	std::cout << "Done." << std::endl;

	size_t nStartTimes = 61;
	arma::vec startTimes = arma::linspace(-2.0, 1.0, nStartTimes); // Make projection times -1:0.05:1

	arma::cube SSIMSproj(nTrials, nDims, nStartTimes, arma::fill::zeros);

	// For each start time, get SSIMS projection
 	std::cout << "Calculate projections: "<< std::endl;
	t1 = high_resolution_clock::now();
	#pragma omp parallel for schedule(dynamic) shared(SSIMSproj)
	for (size_t i = 0; i < nStartTimes; ++i) {
		std::vector<SpikeTrainList_t> projSpkt = getSpikeTrainsFromDataWithWindow(spks, cue_list[2], startTimes[i], winLen);
		arma::mat tmpDmat(nTrials, nBaseSpikeTrains * nEvents);
		for (size_t j = 0; j < nBaseSpikeTrains; ++j) {
			tmpDmat.submat(0, j * nEvents, nTrials - 1, (j + 1) * nEvents - 1) = getFullDistMat(projSpkt[j], basespiketrains[j], myQ);
		}
		SSIMSproj.slice(i) = T1.projectDataIntoTSNESpace(tmpDmat);
	}
	t2 = high_resolution_clock::now();
	std::cout << "Done. Projection of " << nStartTimes << " time steps took: " << std::chrono::duration_cast<std::chrono::milliseconds>( t2 - t1 ).count() << "ms\n";

	size_t i = 0;
	std::cout << "Projection at time " << std::setprecision(3) << startTimes[i] << ":\n" << SSIMSproj.slice(i).rows(0, 3) << " ... " << std::endl << std::endl;

	return 0;
}

int test_data(const char *file)
{
	std::cout << "Reading data from file " << file << std::endl;
	SpkDataList_t cue_list(3);
	SpkData_t trial_group;
	SpkDataList_t spikes;

	high_resolution_clock::time_point t1, t2;

	load_data(file, cue_list, trial_group, spikes);

	t1 = high_resolution_clock::now();
	test_spike_tSNE_move(spikes, cue_list, trial_group);
	t2 = high_resolution_clock::now();
	std::cout << "3 dim tSNE, projection, 10 dim tSNE + classification took " << std::setprecision(0) << std::chrono::duration_cast<std::chrono::milliseconds>( t2 - t1 ).count() << " ms\n";

	t1 = high_resolution_clock::now();
	test_spike_tSNE_combined(spikes, cue_list, trial_group);
	t2 = high_resolution_clock::now();
	std::cout << "5 dim combined tSNE, projection " << std::setprecision(0) << std::chrono::duration_cast<std::chrono::milliseconds>( t2 - t1 ).count() << " ms\n";


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
