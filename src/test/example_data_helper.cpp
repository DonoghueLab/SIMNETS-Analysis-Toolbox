#include <vector>
#include <iostream>
#include <iomanip>
#include <chrono>
using namespace std::chrono;

#include "spikeTrainFun.hpp"
#include "mat.h"

#include "example_data_helper.hpp"

int load_data(const char *file, SpkDataList_t& cue_list, SpkData_t& trial_group, SpkDataList_t& spike_data)
{
	high_resolution_clock::time_point t1, t2;
	t1 = high_resolution_clock::now();
	MATFile *pmat;
	pmat = matOpen(file, "r");
	if (pmat == NULL) {
		std::cout << "Error opening file " << file << std::endl;
		return 1;
	}
	mxArray *spikes = matGetVariable(pmat, "spike_timestamps");
	mxArray *instruction_cue = matGetVariable(pmat, "instruction_cue");
	mxArray *go_cue = matGetVariable(pmat, "go_cue");
	mxArray *start_of_movement = matGetVariable(pmat, "start_of_movement");
	mxArray *movement_direction = matGetVariable(pmat, "movement_direction");

	if (spikes == NULL || instruction_cue == NULL || go_cue == NULL ||
		start_of_movement == NULL || movement_direction == NULL) {
		std::cout << "File " << file << " does not contain the expected variables." << std::endl;
		return 3;
	}

	std::cout << "Spike trains loaded with " << mxGetNumberOfElements(spikes) << " elements." << std::endl;
	std::cout << "Instruction cue times loaded with " << mxGetNumberOfElements(instruction_cue) << " elements." << std::endl;
	std::cout << "Go cue times loaded with " << mxGetNumberOfElements(go_cue) << " elements." << std::endl;
	std::cout << "Start of movement times loaded with " << mxGetNumberOfElements(start_of_movement) << " elements." << std::endl;
	std::cout << "Movement directions loaded with " << mxGetNumberOfElements(movement_direction) << " elements." << std::endl;


	double * data_p = (double*)mxGetData(instruction_cue);
	cue_list[0].insert(cue_list[0].end(), data_p, data_p + mxGetNumberOfElements(instruction_cue));

	data_p = (double*)mxGetData(go_cue);
	cue_list[1].insert(cue_list[1].end(), data_p, data_p + mxGetNumberOfElements(go_cue));

	data_p = (double*)mxGetData(start_of_movement);
	cue_list[2].insert(cue_list[2].end(), data_p, data_p + mxGetNumberOfElements(start_of_movement));

	data_p = (double*)mxGetData(movement_direction);
	trial_group.insert(trial_group.end(), data_p, data_p + mxGetNumberOfElements(movement_direction));

	size_t n_spk = mxGetNumberOfElements(spikes);
	spike_data.resize(n_spk);
	mxArray *tCell;
	for (size_t i = 0; i < n_spk; ++i) {
		tCell = mxGetCell(spikes, i);
		data_p = (double*)mxGetData(tCell);
		spike_data[i].insert(spike_data[i].end(), data_p, data_p + mxGetNumberOfElements(tCell));
	}
	sanitizeSpikeData(spike_data);

	size_t n1 = cue_list[0].size(), n2 = cue_list[1].size(), n3 = cue_list[2].size(), n4 = trial_group.size();
	if (n1!=n2 || n1!=n3 || n1!=n4) {
		std::cout << "event lists should have same length! " << n1 << "; " << n2 << "; " << n3 <<
			"; " << n4 << std::endl;
		return 2;
	}

	mxDestroyArray(spikes);
	mxDestroyArray(instruction_cue);
	mxDestroyArray(go_cue);
	mxDestroyArray(start_of_movement);
	mxDestroyArray(movement_direction);
	t2 = high_resolution_clock::now();
	std::cout << "time elapsed: "<<std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count() / 1e6 << "s" << std::endl;

	return 0;
}

void print_sample_structure(const SpkDataList_t cue_list, const SpkData_t trial_group, const SpkDataList_t spike_data)
{
	std::cout << "Directions and cue times for trials:" << std::endl << std::fixed;
	for(size_t i = 0; i < cue_list[0].size(); ++i) {
		std::cout << std::setprecision(0) << std::setw(5) << trial_group[i] << std::setprecision(3) << std::setw(12) << cue_list[0][i] << std::setw(12) << cue_list[1][i] << std::setw(12) << cue_list[2][i] << std::endl;
	}

	std::cout<< spike_data.size() << " Spike trains, lengths: " << std::endl << std::setprecision(0);
	for (auto c: spike_data) {
		std::cout << std::setw(8) << c.size() << ',';
	}
	std::cout << std::endl;
}
