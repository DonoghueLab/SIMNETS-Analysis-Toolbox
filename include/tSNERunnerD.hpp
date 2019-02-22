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
 * Definition of tSNERunnerD class
 *
 */
#ifndef TSNERUNNERD_HPP_501D143D
#define TSNERUNNERD_HPP_501D143D
#include "tSNERunner.hpp"

/** \brief Direct tSNE transformation

 This class performs tSNE directly on a pair-wise distance matrix of elements.

© Copyright 2016  - Jonas B Zimmermann. All Rights Reserved.

@author Jonas B Zimmermann
@date 2016-08-29
@sa
**/
class TSNERunnerD:public TSNERunner {
public:
	TSNERunnerD(): _pw_dist(defaultMat) {};
	/**	Construct TSNERunnerD with default parameters.
	 *
	 *  @param pw_dist pairwise distance matrix between elements
	 * 	@return TSNERunnerD object.
	**/
	TSNERunnerD(const arma::mat& pw_dist);

	/**	Construct TSNERunnerD with perplexity.
	 *
	 *  @param pw_dist pairwise distance matrix between elements
	 * 	@param perplexity defines the notion of neighbourhood for tSNE.
	 * 	@return TSNERunnerD object.
	**/
	TSNERunnerD (const arma::mat& pw_dist, double perplexity);

	/**	Construct TSNERunnerD with number of dimensions and perplexity.
	 *
	 *  @param pw_dist pairwise distance matrix between elements
	 * 	@param outDimensions number of dimensions of tSNE space.
	 * 	@param perplexity defines the notion of neighbourhood for tSNE.
	 * 	@return TSNERunnerD object.
	**/
	TSNERunnerD (const arma::mat& pw_dist, size_t outDimensions, double perplexity);

	/**	Construct TSNERunnerD with number of dimensions and perplexity.
	 *
	 *  @param pw_dist pairwise distance matrix between elements
	 * 	@param outDimensions number of dimensions of tSNE space.
	 * 	@return TSNERunnerD object.
	**/
	TSNERunnerD (const arma::mat& pw_dist, size_t outDimensions);

	/**	Construct TSNERunnerD with default parameters.
	 *
	 *  @param pw_dist pairwise distance matrix between elements
	 * 	@param params: a TSNERunnerAlgoParams struct containing parameters used in the tSNE algorithm.
	 * 	@return TSNERunnerD object.
	**/
	TSNERunnerD(const arma::mat& pw_dist, TSNERunnerAlgoParams params);

	/**	Construct TSNERunnerD with perplexity.
	 *
	 *  @param pw_dist pairwise distance matrix between elements
	 * 	@param params: a TSNERunnerAlgoParams struct containing parameters used in the tSNE algorithm.
	 * 	@param perplexity defines the notion of neighbourhood for tSNE.
	 * 	@return TSNERunnerD object.
	**/
	TSNERunnerD (const arma::mat& pw_dist, double perplexity, TSNERunnerAlgoParams params);

	/**	Construct TSNERunnerD with number of dimensions and perplexity.
	 *
	 *  @param pw_dist pairwise distance matrix between elements
	 * 	@param params: a TSNERunnerAlgoParams struct containing parameters used in the tSNE algorithm.
	 * 	@param outDimensions number of dimensions of tSNE space.
	 * 	@param perplexity defines the notion of neighbourhood for tSNE.
	 * 	@return TSNERunnerD object.
	**/
	TSNERunnerD (const arma::mat& pw_dist, size_t outDimensions, double perplexity, TSNERunnerAlgoParams params);

	/**	Construct TSNERunnerD with number of dimensions and perplexity.
	 *
	 *  @param pw_dist pairwise distance matrix between elements
	 * 	@param params: a TSNERunnerAlgoParams struct containing parameters used in the tSNE algorithm.
	 * 	@param outDimensions number of dimensions of tSNE space.
	 * 	@return TSNERunnerD object.
	**/
	TSNERunnerD (const arma::mat& pw_dist, size_t outDimensions, TSNERunnerAlgoParams params);

	virtual ~TSNERunnerD() {};

	/**	@brief Run tSNE on data
	 *
	 *	Perform tSNE dimensionality reduction on data.
	 *	@return Kullback-Leibler divergence
	**/
	double run (void) override;

protected:
	const arma::mat& _pw_dist; ///< Stores pairwise distance data
	static arma::mat defaultMat;
};

#endif /* end of include guard: TSNERUNNERD_HPP_501D143D */
