/**
 * @file spikeTrainFun.hpp
 * @author  Jonas Zimmermann <jonas_zimmermann@brown.edu>
 * @version 0.5
 * @brief Defines class to handle spike train data, as well as some helper functions to calculate spike train distances for lists of spike trains.
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
 * Definition of SpikeTrain class, which stores spike times and offers basic
 * operations, such as distances. Further, this file includes helper functions
 * calculate distance matrices between spike trains.
 */
#ifndef SPIKETRAINFUN_HPP_C855D429
#define SPIKETRAINFUN_HPP_C855D429

#include <array>
#include <vector>
#include <armadillo>
#include <memory>

class SpikeTrain;

/** \brief Data type for spike times.
 * Could also be float, or an integer type (for ticks)
**/
typedef double SpkDatum_t;

/** \brief Container for multiple spikes.
 * Stores raw data for spike trains
**/
typedef std::vector<SpkDatum_t> SpkData_t;

/** \brief Container for multiple raw spike trains.
 * To store several raw spike trains, for example from several units
**/
typedef std::vector<SpkData_t> SpkDataList_t;

/** \brief Container for shared_ptrs to multiple raw spike trains
 * 	To store pointers to several raw spike trains
**/
typedef std::vector<std::shared_ptr<SpkData_t>> SpkDataPtrList_t;

/** \brief Container for multiple SpikeTrain objects.
**/
typedef std::vector<SpikeTrain> SpikeTrainList_t;

/** \brief Container for spike train distances.
**/
typedef std::vector<double> Dist_t;

/** \brief Class to store spike trains.
 *
 * A spike train stores a reference to the actual spike train
 * data and a window that determines which spikes to be considered for
 * calculations such as the Victor&Purpura spike train metric. This
 * implementation ensures that only a minimum of data is stored, making
 * calculations based on spike trains derived from the same underlying data
 * fast and efficient.
**/
class SpikeTrain
{
public:
/**
 * Constructor that creates an empty spike train.
**/
	SpikeTrain (void);

/**
 * Constructor that creates a spike train with spike data and an all-
 * encompassing window.
 *
 * @param t \ref SpkData_t containing spike times
**/
	SpikeTrain (const SpkData_t& t);

/**
 * Constructor that creates a spike train with data and a window.
 *
 * @param t \ref SpkData_t containing spike times
 * @param win_t array of 2 \ref SpkDatum_t, marking beginning and end of window
**/
	SpikeTrain (const SpkData_t& t, const std::array<SpkDatum_t, 2> win_t);

	/**
	* Constructor that creates a spike train with data and a window.
	*
	* @param t \ref SpkData_t containing spike times
	* @param win_t array of 2 \ref SpkDatum_t, marking beginning and end of window
	* @param offset \ref SpkDatum_t, offset of spike train times relative to window
	**/
	SpikeTrain (const SpkData_t& t, const std::array<SpkDatum_t, 2> win_t, const SpkDatum_t offset);

/**
 * Constructor that creates a spike train with a pointer to spike data and an all-
 * encompassing window.
 *
 * @param t_ptr std::shared_ptr to \ref SpkData_t containing spike times
**/
	SpikeTrain (const std::shared_ptr<SpkData_t> t_ptr);

/**
 * Constructor that creates a spike train with a pointer to data and a window.
 *
 * @param t_ptr std::shared_ptr to \ref SpkData_t containing spike times
 * @param win_t array of 2 \ref SpkDatum_t, marking beginning and end of window
**/
	SpikeTrain (const std::shared_ptr<SpkData_t> t_ptr, const std::array<SpkDatum_t, 2> win_t);

   /**
    * Constructor that creates a spike train with a pointer to data and a window.
    *
    * @param t_ptr std::shared_ptr to \ref SpkData_t containing spike times
    * @param win_t array of 2 \ref SpkDatum_t, marking beginning and end of window
	* @param offset \ref SpkDatum_t, offset of spike train times relative to window
   **/
   	SpikeTrain (const std::shared_ptr<SpkData_t> t_ptr, const std::array<SpkDatum_t, 2> win_t, const SpkDatum_t offset);

	/**
	* Constructor that creates a spike train with a pointer to data and a window.
	*
	* @param t_ptr std::shared_ptr to \ref SpkData_t containing spike times
	* @param winStart \ref SpkDatum_t denoting beginning of window
	* @param winLength \ref SpkDatum_t denoting length of window
	**/
	SpikeTrain (const std::shared_ptr<SpkData_t> t_ptr, const SpkDatum_t winStart, const SpkDatum_t winLength);
	/**
	* Constructor that creates a spike train with a pointer to data and a window.
	*
	* @param t_ptr std::shared_ptr to \ref SpkData_t containing spike times
	* @param winStart \ref SpkDatum_t denoting beginning of window
	* @param winLength \ref SpkDatum_t denoting length of window
	* @param offset \ref SpkDatum_t, offset of spike train times relative to window
	**/
	SpikeTrain (const std::shared_ptr<SpkData_t> t_ptr, const SpkDatum_t winStart, const SpkDatum_t winLength, const SpkDatum_t offset);
/**
 * Constructor that creates a spike train from another one.
 *
 * @param other SpikeTrain object
**/
	SpikeTrain (const SpikeTrain& other);

/**
 * Returns use_count of std::shared_ptr \ref t_ .
**/
	long int spkUseCount() const {return t_.use_count();}

/**
 * Returns spikes included in window as raw spike data.
 * @return raw spike data
**/
	SpkData_t spikesInWindow(void) const;

/**
 * Returns number of spikes included in window.
**/
	SpkData_t::size_type length(void) const {return std::distance(win_i_[0], win_i_[1]);}

/**
 * Returns number of spikes included in window.
**/
	SpkData_t::size_type size(void) const {return length();}

/**
 * Returns std::weak_ptr to raw underlying spike data.
**/
	std::weak_ptr<SpkData_t> getSpikes() {return t_;}

/**
 * Print SpikeTrain, including window.
**/
	void print (void) const;

/**
 * Print SpikeTrain, including window.
 * @param os Output stream
**/
	void print (std::ostream& os) const;

/**
 * Calculate Victor & Purpura distance between spike trains.
 * @param otherSt SpikeTrain to compare to
 * @param q Cost for moving a spike (temporal resolution)
 * @return VP distance
**/
	double distVPTo(const SpikeTrain& otherSt, const double q) const;

/**
 * Set a new window for spike train.
 * @param win array of 2 \ref SpkDatum_t, marking beginning and end of window
**/
	void setWin(const std::array<SpkDatum_t, 2>& win);

/**
 * Set a new window for spike train.
 * @param winStart \ref SpkDatum_t denoting beginning of window
 * @param winLength \ref SpkDatum_t denoting length of window
**/
	void setWin(const SpkDatum_t winStart, const SpkDatum_t winLength);

/**
 * Apply window (stored as beginning and end times). Finds first and last spike
 * in window (if any) and stores them as iterators in \ref win_i_.
**/
	void applyWin(void);

/**
 * Determines whether a time falls into window or not.
 * @param s spike time to test
 * @return true if s falls into window, false otherwise
**/
	inline bool sInWin(SpkDatum_t s) const { return ((win_t_[0] <= s) && (s < win_t_[1]));}
private:
	std::shared_ptr<SpkData_t> t_; /**< stores a reference to the actual spike train data **/
	std::array<SpkDatum_t, 2> win_t_; /**< stores the window imposed on the spike data **/
	std::array<SpkData_t::iterator, 2> win_i_; /**< stores the window imposed on the spike data as iterators on \ref t_ **/
	SpkDatum_t _offset = 0; /**< offset applied to spike times when reported relative to window **/
};

/**
 * Overload << operator to allow printing of SpikeTrain objects
 * @param os std::ostream Output stream
 * @param st SpikeTrain object
 * @return std::ostream
**/
std::ostream& operator<<(std::ostream& os, const SpikeTrain& st);

/**
 * Function to generate a number SpikeTrain objects by applying a window
 * to a list of raw spike train data.
 * @param tl list of raw spike trains
 * @param win array of 2 \ref SpkDatum_t, marking beginning and end of window
 * @return \ref SpikeTrainList_t vector of SpikeTrain objects; length equal to length of tl
**/
SpikeTrainList_t getSpikeTrainsFromDataWithWindow(const SpkDataList_t& tl, const std::array<SpkDatum_t, 2> win);

/**
 * Function to generate a number SpikeTrain objects by applying a window
 * to a list of raw spike train data
 * @param tl list of raw spike trains
 * @param winStart \ref SpkDatum_t denoting beginning of window
 * @param winLength \ref SpkDatum_t denoting length of window
 * @return \ref SpikeTrainList_t vector of SpikeTrain objects; length equal to length of tl
**/
SpikeTrainList_t getSpikeTrainsFromDataWithWindow(const SpkDataList_t& tl, const SpkDatum_t winStart, const SpkDatum_t winLength);

/**
 * Function to generate a number SpikeTrain objects by applying a window
 * to raw spike train data, around a number of time points.
 * @param t raw spike train
 * @param aTimes time points around which to place window
 * @param win array of 2 \ref SpkDatum_t, marking beginning and end of window
 * @return \ref SpikeTrainList_t vector of SpikeTrain objects; length equal to length of aTimes
**/
SpikeTrainList_t getSpikeTrainsFromDataWithWindow(const SpkData_t& t, const SpkData_t& aTimes, const std::array<SpkDatum_t, 2> win);

/**
 * Function to generate a number SpikeTrain objects by applying a window
 * to raw spike train data, around a number of time points.
 * @param t raw spike train
 * @param aTimes time points around which to place window
 * @param winStart \ref SpkDatum_t denoting beginning of window
 * @param winLength \ref SpkDatum_t denoting length of window
 * @return \ref SpikeTrainList_t vector of SpikeTrain objects; length equal to length of aTimes
**/
SpikeTrainList_t getSpikeTrainsFromDataWithWindow(const SpkData_t& t, const SpkData_t& aTimes, const SpkDatum_t winStart, const SpkDatum_t winLength);

/**
 * Function to generate a number SpikeTrain objects by applying a window
 * to raw spike train data, around a number of time points.
 * @param tl list of raw spike trains
 * @param aTimes time points around which to place window
 * @param winStart \ref SpkDatum_t denoting beginning of window
 * @param winLength \ref SpkDatum_t denoting length of window
 * @return \ref SpikeTrainList_t vector of SpikeTrain objects; length equal to length of aTimes
**/
std::vector<SpikeTrainList_t> getSpikeTrainsFromDataWithWindow(const SpkDataList_t& tl, const SpkData_t& aTimes, const SpkDatum_t winStart, const SpkDatum_t winLength);


SpikeTrainList_t getSpikeTrainsFromDataWithWindow(const std::shared_ptr<SpkData_t>& t, const SpkData_t& aTimes, const std::array<SpkDatum_t, 2> win);
SpikeTrainList_t getSpikeTrainsFromDataWithWindow(const std::shared_ptr<SpkData_t>& t, const SpkData_t& aTimes, const SpkDatum_t winStart, const SpkDatum_t winLength);
SpikeTrainList_t getSpikeTrainsFromDataWithWindow(const SpkDataPtrList_t& tl, const std::array<SpkDatum_t, 2> win);
SpikeTrainList_t getSpikeTrainsFromDataWithWindow(const SpkDataPtrList_t& tl, const SpkDatum_t winStart, const SpkDatum_t winLength);
std::vector<SpikeTrainList_t> getSpikeTrainsFromDataWithWindow(const SpkDataPtrList_t& tl, const SpkData_t& aTimes, const SpkDatum_t winStart, const SpkDatum_t winLength);

/**
 * Calculate spike counts in event windows
 * @param spikeTrains \ref SpkDataList_t vector of $n$ SpkData_t objects (i.e. vector<vector<double>>)
 * @param events \ref SpkDataList_t vector of $m$ event lists relative to which spikes will be counted
 * @param offsets \ref SpkData_t vector of $m$ offset times (relative to times given in events)
 * @param lengths \ref SpkData_t vector of $m$ durations for time windows (after event[j] + offsets[j])
 * @return vector of vectors of vectors of size_t: $n$ spike trains with $m$ vectors each representing event list, each of those having $k_j$ elements ($k_j$ is length of $j$-th event vector)
**/
std::vector<std::vector<std::vector<size_t>>> getSpikeDataEventWindowCounts(const SpkDataList_t &spikeTrains, const SpkDataList_t &events, const SpkData_t &offsets, const SpkData_t &lengths);

/**
 * Return spike times in event windows
 * @param spikeTrains \ref SpkDataList_t vector of $n$ SpkData_t objects (i.e. vector<vector<double>>)
 * @param events \ref SpkDataList_t vector of $m$ event lists relative to which spikes will be counted
 * @param offsets \ref SpkData_t vector of $m$ offset times (relative to times given in events)
 * @param lengths \ref SpkData_t vector of $m$ durations for time windows (after event[j] + offsets[j])
 * @return vector of vectors of vectors of double: $n$ spike trains with $m$ vectors each representing event list, each of those having $k_j$ elements ($k_j$ is length of $j$-th event vector)
**/
std::vector<std::vector<SpkDataList_t>> getSpikeDataInEventWindows(const SpkDataList_t &spikeTrains, const SpkDataList_t &events, const SpkData_t &offsets, const SpkData_t &lengths);

/**
 * Calculate pairwise spike train distances.
 * @param spikeTrains \ref SpikeTrainList_t vector of n SpikeTrain objects
 * @param q Cost for moving a spike (temporal resolution)
 * @return \ref Dist_t vector of distances, contains \f$n(n - 1)/2\f$ elements:
 * according to scheme (1, 2), (1, 3), ... ,(1, n), (2, 3), ..., (n - 1, n)
**/
Dist_t getDistMat(const SpikeTrainList_t& spikeTrains, const double q);

/**
 * Calculate pairwise spike train distances. Pass iterator pointing to some Dist_t
 * object which should be large enough to hold all \f$n (n - 1) /2\f$ values,
 * where \f$n\f$ is the number of spike trains
 * @param spikeTrains \ref SpikeTrainList_t vector of n SpikeTrain objects
 * @param q Cost for moving a spike (temporal resolution)
 * @param dBegin iterator pointing to start where distance data will be written to
 * @param dEnd iterator pointing to end of Dist_t.
 * @return \ref Dist_t iterator pointing to the next element in Dist_t that would be written to.
**/
Dist_t::iterator getDistMat(const SpikeTrainList_t& spikeTrains, const double q,
	Dist_t::iterator dBegin, const  Dist_t::const_iterator dEnd);

arma::mat getFullDistMat(const SpikeTrainList_t& spikeTrains, const double q);

/**
 * Calculate pairwise spike train distances between two sets of SpikeTrain objects.
 * @param spikeTrains1 \ref SpikeTrainList_t vector of n SpikeTrain objects
 * @param spikeTrains2 \ref SpikeTrainList_t vector of m SpikeTrain objects
 * @param q Cost for moving a spike (temporal resolution)
 * @return \ref Dist_t vector of distances, contains \f$n\cdot m\f$ elements:
 * according to scheme (1, 1), (2, 1), ... ,(n, 1), (1, 2), ..., (n, m)
**/
Dist_t getDistMat(const SpikeTrainList_t& spikeTrains1, const SpikeTrainList_t& spikeTrains2, const double q);

/**
 * Calculate pairwise spike train distances between two sets of SpikeTrain objects.
 * Pass iterator pointing to some Dist_t object which should be large enough to
 * hold all \f$n m\f$ values, where \f$n\f$ is the number in the first set and
 * \f$m\f$ is the number in the second set
 * @param spikeTrains1 \ref SpikeTrainList_t vector of n SpikeTrain objects
 * @param spikeTrains2 \ref SpikeTrainList_t vector of m SpikeTrain objects
 * @param q Cost for moving a spike (temporal resolution)
 * @param dBegin iterator pointing to start where distance data will be written to
 * @param dEnd iterator pointing to end of Dist_t.
 * @return \ref Dist_t iterator pointing to the next element in Dist_t that would be written to.
**/
Dist_t::iterator getDistMat(const SpikeTrainList_t& spikeTrains1, const SpikeTrainList_t& spikeTrains2,
	const double q, Dist_t::iterator dBegin, const Dist_t::const_iterator dEnd);

arma::mat getFullDistMat(const SpikeTrainList_t& spikeTrains1, const SpikeTrainList_t& spikeTrains2, const double q);

/**
 * Calculate spike train distances between neurons, at various time points.
 * Pass iterator pointing to some \ref Dist_t object which should be large enough to
 * hold all \f$n\cdot m\cdot (m - 1) /2\f$ values, where \f$m\f$ is the number
 * of \ref SpikeTrain objects in tl, and \f$n\f$ is the number of time points in aTimes.
 * Pass another set of iterators for \ref SpkDataList_t to hold \f$n\cdot m\f$ base \ref SpikeTrain objects.
 * @param tl list of raw spike trains
 * @param aTimes time points around which to place window
 * @param winStart \ref SpkDatum_t denoting beginning of window
 * @param winLength \ref SpkDatum_t denoting length of window
 * @param q Cost for moving a spike (temporal resolution)
 * @param dBegin iterator pointing to start where distance data will be written to
 * @param dEnd iterator pointing to end of \ref Dist_t.
 * @param sBegin iterator pointing to start where base spiketrains will be written to
 * @param sEnd iterator pointing to end of \ref SpkDataList_t.
**/
void getSSIMDMatBetweenNeurons(const SpkDataList_t& tl, const SpkData_t& aTimes,
	const SpkDatum_t winStart, const SpkDatum_t winLength, const double q,
	Dist_t::iterator dBegin, const Dist_t::const_iterator dEnd,
	SpkDataList_t::iterator sBegin, const SpkDataList_t::const_iterator sEnd);

arma::mat getSSIMDMatBetweenNeurons(const SpkDataList_t& tl, const SpkData_t& aTimes,
	const SpkDatum_t winStart, const SpkDatum_t winLength, const double q,
	SpkDataList_t::iterator sBegin, const SpkDataList_t::const_iterator sEnd);

arma::mat getSSIMDMatBetweenNeurons(const SpkDataPtrList_t& tl, const SpkData_t& aTimes,
	const SpkDatum_t winStart, const SpkDatum_t winLength, const double q,
	SpkDataList_t::iterator sBegin, const SpkDataList_t::const_iterator sEnd);


/**
 * Calculate spike train distances between time points, for several neurons.
 * Pass iterator pointing to some \ref Dist_t object which should be large enough to
 * hold all \f$n\cdot m\cdot (m - 1) /2\f$ values, where \f$n\f$ is the number
 * of \ref SpikeTrain objects in tl, and \f$m\f$ is the number of time points in aTimes.
 * Pass another set of iterators for \ref SpkDataList_t to hold \f$n\cdot m\f$ base \ref SpikeTrain objects.
 * @param tl list of raw spike trains
 * @param aTimes time points around which to place window
 * @param winStart \ref SpkDatum_t denoting beginning of window
 * @param winLength \ref SpkDatum_t denoting length of window
 * @param q Cost for moving a spike (temporal resolution)
 * @param dBegin iterator pointing to start where distance data will be written to
 * @param dEnd iterator pointing to end of \ref Dist_t.
 * @param sBegin iterator pointing to start where base spiketrains will be written to
 * @param sEnd iterator pointing to end of \ref SpkDataList_t.
**/
void getSSIMDMatBetweenTimePoints(const SpkDataList_t& tl, const SpkData_t& timePoints,
	const SpkDatum_t winStart, const SpkDatum_t winLength, const double q,
	Dist_t::iterator dBegin, const Dist_t::const_iterator dEnd,
	SpkDataList_t::iterator sBegin, const SpkDataList_t::const_iterator sEnd);

/**
* Calculate spike train distances between time points, for several neurons.
* Pass iterator pointing to some \ref Dist_t object which should be large enough to
* hold all \f$n\cdot m\cdot (m - 1) /2\f$ values, where \f$n\f$ is the number
* of \ref SpikeTrain objects in tl, and \f$m\f$ is the number of time points in aTimes.
* Pass another set of iterators for \ref SpkDataList_t to hold \f$n\cdot m\f$ base \ref SpikeTrain objects.
* @param tl list of raw spike trains
* @param aTimes time points around which to place window
* @param winStart \ref SpkDatum_t denoting beginning of window
* @param winLength \ref SpkDatum_t denoting length of window
* @param q Cost for moving a spike (temporal resolution)
* @param sBegin iterator pointing to start where base spiketrains will be written to
* @param sEnd iterator pointing to end of \ref SpkDataList_t.
**/
arma::mat getSSIMDMatBetweenTimePoints(const SpkDataList_t& tl, const SpkData_t& timePoints,
const SpkDatum_t winStart, const SpkDatum_t winLength, const double q,
SpkDataList_t::iterator sBegin, const SpkDataList_t::const_iterator sEnd);

arma::mat getSSIMDMatBetweenTimePoints(const SpkDataPtrList_t& tl, const SpkData_t& timePoints,
const SpkDatum_t winStart, const SpkDatum_t winLength, const double q,
SpkDataList_t::iterator sBegin, const SpkDataList_t::const_iterator sEnd);


/**
 * Copy pairwise distance data to a square matrix.
 * @param dm pairwise distance vector of length \f$n\f$
 * @param dest pointer to array of double, which has to be large enough to contain all
**/
size_t copyPDistToSquare(const Dist_t &dm, double *dest);

arma::mat copyNTimesPDistToMat(const Dist_t &dm, const size_t n_reps, const size_t n_pw);

/**
 * Sanitize raw spike data, i.e. sort spike times
 * @param tl list of raw spike trains
**/
void sanitizeSpikeData(SpkDataList_t& tl);

/**
 * Sanitize raw spike data, i.e. sort spike times
 * @param t raw spike train
**/
void sanitizeSpikeData(SpkData_t& t);

/**
 * Sanitize raw spike data, i.e. sort spike times
 * @param tl vector of shared_ptrs to raw spike trains
**/
void sanitizeSpikeData(SpkDataPtrList_t& tl);

/**
 * Sanitize raw spike data, i.e. sort spike times
 * @param t shared_ptr to raw spike train
**/
void sanitizeSpikeData(std::shared_ptr<SpkData_t> t);


/** Load spike trains from binary data
 *
**/
SpkDataPtrList_t loadSpikeTrains(const char* fn);

#endif /* end of include guard: SPIKETRAINFUN_HPP_C855D429 */
