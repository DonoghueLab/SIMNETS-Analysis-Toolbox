#ifndef EXAMPLE_DATA_HELPER_HPP_BA0D6C8E
#define EXAMPLE_DATA_HELPER_HPP_BA0D6C8E

//#include "spikeTrainFun.hpp"

int load_data(const char *file, SpkDataList_t& cue_list, SpkData_t& trial_group, SpkDataList_t& spike_data);
void print_sample_structure(const SpkDataList_t cue_list, const SpkData_t trial_group, const SpkDataList_t spike_data);

#endif /* end of include guard: EXAMPLE_DATA_HELPER_HPP_BA0D6C8E */
