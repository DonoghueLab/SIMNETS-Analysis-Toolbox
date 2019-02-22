#ifndef SSIMSHELPER_HPP_33B169F2
#define SSIMSHELPER_HPP_33B169F2
/**
 * @file ssimsHelper.hpp
 * @author  Jonas Zimmermann <jonas_zimmermann@brown.edu>
 * @version 0.5
 * @brief Defines helper functions for SSIMS transformations.
 *
 * @section LICENSE
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
 * @section DESCRIPTION
 *
 * Use functions defined here in order to facilitate SSIMS analyses.
**/

/** \brief Run sweep over various SSIMS parameters
 *
 * 	This function runs SSIMS for a given collection of neurons and events and returns the low-dimensional SSIMS
 * 	spaces at every parameter combination
 *
 * \author Jonas B Zimmermann
 * \date 2015-11-24
 * \param  tl: vector of spike trains
 * \param  evTimes: time points of events relative to which windows will be anchored
 * \param  eventOffsets: first parameter of sweep; offset to add to evTimes for windows. Minor index for output
 * \param  winLens: 2nd parameter of sweep; length of spike train windows
 * \param  qs: 3rd parameter of sweep; cost parameter for spike train metric
 * \param  nDims: 4th parameter of sweep; number of output dimensions
 * \param  perplexities: 5th parameter of sweep; perplexity values to use for tSNE. Major index for output
 * \param  direct_tSNE if true, performs direct tSNE (i.e. do not consider pw distance matrix as feature vector) if tl is a single spike train
 * \return vector of SSIMS matrices, indexed by sweep parameters: n_eo * (n_wl * (n_q * i_p + i_q) + i_wl) + i_eo
 * \sa
**/
std::vector<arma::mat> doSsimsSweep(const SpkDataList_t &tl, const SpkData_t &evTimes, const SpkData_t &eventOffsets, const SpkData_t &winLens, const SpkData_t &qs, const std::vector<size_t> &nDims, const SpkData_t &perplexities, bool direct_tSNE);


#endif /* end of include guard: SSIMSHELPER_HPP_33B169F2 */
