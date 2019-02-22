#include <armadillo>
#include <vector>
#include <iterator>
#include "tSNERunnerD.hpp"
#include "tSNERunnerV.hpp"
#include "spikeTrainFun.hpp"
#include "ssimsHelper.hpp"

std::vector<arma::mat> doSsimsSweep(const SpkDataList_t &tl, const SpkData_t &evTimes, const SpkData_t &eventOffsets, const SpkData_t &winLens, const SpkData_t &qs, const std::vector<size_t> &nDims, const SpkData_t &perplexities, bool direct_tSNE)
{
	size_t nEventOffsets = eventOffsets.size(), nWinLens = winLens.size(), nQs = qs.size();
	size_t nNDims = nDims.size(), nPerplexities = perplexities.size();
	size_t nSpikeTrains = tl.size();

	std::vector<arma::mat> outSSIMS(nEventOffsets * nWinLens * nQs * nNDims * nPerplexities);
	SpkDataList_t baseSpkTr;

	#pragma omp parallel for shared(outSSIMS)
	for (size_t i_eo = 0; i_eo < nEventOffsets; ++i_eo) {
		for (size_t i_wl = 0; i_wl < nWinLens; ++i_wl) {
			for (size_t i_q = 0; i_q < nQs; ++i_q) {
				for (size_t i_d = 0; i_d < nNDims; ++i_d) {
					for (size_t i_p = 0; i_p < nPerplexities; ++i_p) {
						size_t ind = nEventOffsets * (nWinLens * ( nQs * ( nNDims * i_p + i_d) + i_q) + i_wl) + i_eo;
						SpkDatum_t evOffset = eventOffsets[i_eo], winLen = winLens[i_wl], q = qs[i_q];
						SpkDatum_t perplexity = perplexities[i_p];
						size_t nDim = nDims[i_d];
						arma::mat aBaseData = getSSIMDMatBetweenTimePoints( tl, evTimes,
						evOffset, winLen, q, baseSpkTr.begin(), baseSpkTr.end());
						TSNERunner *tSNE;
						if (direct_tSNE && (nSpikeTrains == 1)) {
							tSNE = new TSNERunnerD(aBaseData, nDim, perplexity);
						}
						else {
							tSNE = new TSNERunnerV(aBaseData, nDim, perplexity);
						}
						tSNE->run();
						outSSIMS[ind] = tSNE->getTSNEData();
						delete tSNE;
					}
				}
			}
		}
	}
	return outSSIMS;
}
