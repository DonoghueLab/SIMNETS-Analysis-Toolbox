/**
 * @file tSNERunner.hpp
 * @author  Jonas Zimmermann <jonas_zimmermann@brown.edu>
 * @version 0.01
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
 * Definition of tSNERunner class.
 *
 */
#ifndef TSNERUNNER_HPP_FC374924
#define TSNERUNNER_HPP_FC374924
#include <armadillo>

enum NormMethod {	 	///< Defines normalization method by which data will be processed before applying tSNE
	NMNoNorm,				///< No normalization will be applied
	NMByFeature, 			///< Normalization will be applied by feature, i.e. columns of the data matrix
	NMGlobal 				///< A global normalization will be applied
};

/** Defines parameters for the tSNE algorithm.
 *
 * For details regarding these parameters, consult the [original
 * t-SNE paper][1]. The default values are sensible and may not
 * need to be tuned much.
 *
 * [1]: http://www.jmlr.org/papers/v9/vandermaaten08a.html "Van der Maaten, Laurens J P and Geoffrey E Hinton (Nov. 2008). “Visualizing High-Dimensional Data Using t-SNE”. In: Journal of Machine Learning Research 9, pp. 2579–2605."
**/
struct TSNERunnerAlgoParams
{
	TSNERunnerAlgoParams ();
	size_t momSwitchIter;
	size_t stopLyingIter;
	size_t maxIter;
	double epsilon;
	double minGain;
	double momentum;
	double finalMomentum;
};

/** \brief  abstract base class for running t-SNE.

    This is the base class for all tSNE operations. It defines the interface
    for all implementations.

© Copyright 2016  - Jonas B Zimmermann. All Rights Reserved.

@author Jonas B Zimmermann
@date 2016-08-29
@sa
**/
class TSNERunner {
public:
	enum State {
		Uninitialized,
		JointPCalculated,
		Ready
	};


	/**	Construct TSNERunner with default parameters.
	 *
	 * 	@return TSNERunner object.
	**/
	TSNERunner();

	/**	Construct TSNERunner with perplexity.
	 *
	 * 	@param perplexity defines the notion of neighbourhood for tSNE.
	 * 	@return TSNERunner object.
	**/
	TSNERunner (double perplexity);

	/**	Construct TSNERunner with number of dimensions and perplexity.
	 *
	 * 	@param outDimensions number of dimensions of tSNE space.
	 * 	@param perplexity defines the notion of neighbourhood for tSNE.
	 * 	@return TSNERunner object.
	**/
	TSNERunner (size_t outDimensions, double perplexity);

	/**	Construct TSNERunner with number of dimensions
	 *
	 * 	@param outDimensions number of dimensions of tSNE space.
	 * 	@return TSNERunner object.
	**/
	TSNERunner (size_t outDimensions);

	/**	Construct TSNERunner with default parameters.
	 *
	 * 	@param params: a TSNERunnerAlgoParams struct containing parameters used in the tSNE algorithm.
	 * 	@return TSNERunner object.
	**/
	TSNERunner(TSNERunnerAlgoParams params);

	/**	Construct TSNERunner with perplexity.
	 *
	 * 	@param perplexity defines the notion of neighbourhood for tSNE.
	 * 	@param params: a TSNERunnerAlgoParams struct containing parameters used in the tSNE algorithm.
	 * 	@return TSNERunner object.
	**/
	TSNERunner (double perplexity, TSNERunnerAlgoParams params);

	/**	Construct TSNERunner with number of dimensions and perplexity.
	 *
	 * 	@param outDimensions number of dimensions of tSNE space.
	 * 	@param perplexity defines the notion of neighbourhood for tSNE.
	 * 	@param params: a TSNERunnerAlgoParams struct containing parameters used in the tSNE algorithm.
	 * 	@return TSNERunner object.
	**/
	TSNERunner (size_t outDimensions, double perplexity, TSNERunnerAlgoParams params);

	/**	Construct TSNERunner with number of dimensions
	 *
	 * 	@param outDimensions number of dimensions of tSNE space.
	 * 	@param params: a TSNERunnerAlgoParams struct containing parameters used in the tSNE algorithm.
	 * 	@return TSNERunner object.
	**/
	TSNERunner (size_t outDimensions, TSNERunnerAlgoParams params);

	virtual ~TSNERunner() {};

	/// set Perplexity parameter for tSNE
	void setPerplexity (double perplexity);

	/// get current Perplexity parameter
	double getPerplexity (void) const;

	/// Set normalization method during preprocessing
	void setNormalizationMethod(NormMethod preProcNorm);

	/// Get normalization method during preprocessing
	NormMethod getNormalizationMethod (void) const;

	/// Set initial solution to supplied matrix, which should have the same size as the output data (i.e. be a _baseData.n_rows x _outDimensions matrix)
	virtual void setInitialSolution(const arma::mat& initialSolution);


	/**	@brief Run tSNE on data
	 *
	 *	Perform tSNE dimensionality reduction on data.
	 *	@return Kullback-Leibler divergence
	**/
	virtual double run (void) = 0;

	/**	Get tSNE-transformed data
	 *
	 * 	After tSNE has been run, use this function to return the transformed data.
	 * 	@return \f$n\times m\f$-matrix, where \f$n\f$ equals number of rows of _baseData and \f$m\f$ equals _outDimensions.
	**/
	const arma::mat& getTSNEData (void) const;

protected:
	/// Delete any solutions found, reset State
	virtual void reinitialize (void);

	TSNERunner::State _state = Uninitialized;
	NormMethod _preProcNorm = NMGlobal; ///< normalization method by which data will be processed before applying tSNE
	TSNERunnerAlgoParams _params; ///< keeps parameters for the tSNE algorithm
	arma::mat _ydata; ///< Stores tSNE-transformed data
	double _perplexity = 30; ///< Perplexity setting of the tSNE algorithm, roughly translating to number of nearest neighbours

	double _klDiv = 0;
	size_t _outDimensions = 2; ///< Dimensions of the final tSNE embedding

};

arma::mat computeJointProbability(const arma::mat& data, NormMethod nm, double perplexity);

/**
 * Identifies appropriate sigma's to get kk NNs up to some tolerance
 * Adopted directly from van der Maaten's MATLAB implementation
**/

arma::mat d2p(const arma::mat& D, const double u=15, const double tol=1e-4);

std::pair<double, arma::mat> Hbeta(const arma::mat& D, const double beta);

#endif /* end of include guard: TSNERUNNER_HPP_FC374924 */
