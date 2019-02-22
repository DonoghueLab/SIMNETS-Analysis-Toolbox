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
 * Definition of tSNERunnerV class
 *
 */
#ifndef TSNERUNNERV_HPP_938F72DD
#define TSNERUNNERV_HPP_938F72DD
#include "tSNERunner.hpp"

/** \brief Runs tSNE after PCA

 This class allows to run tSNE on data which has more dimensions than observations.
 To this end, a principal component analysis is performed first, reducing the
 dimensionality of the data to the number of observations, then performs tSNE.
 This algorithm then also provides a means to project new data into the tSNE
 space, using the method @ref projectDataIntoTSNESpace.

© Copyright 2016  - Jonas B Zimmermann. All Rights Reserved.

@author Jonas B Zimmermann
@date 2016-08-29
@sa
**/
class TSNERunnerV:public TSNERunner
{
public:
	/**	Construct TSNERunnerV with matrix.
	 *
	 * 	@param baseData is a matrix, observations in rows and variables in columns.
	 * 	@return TSNERunner object.
	**/
	TSNERunnerV (const arma::mat& baseData);

	/**	Construct TSNERunnerV with matrix and perplexity.
	 *
	 * 	@param baseData is a matrix, observations in rows and variables in columns.
	 * 	@param perplexity defines the notion of neighbourhood for tSNE.
	 * 	@return TSNERunner object.
	**/
	TSNERunnerV (const arma::mat& baseData, double perplexity);

	/**	Construct TSNERunnerV with matrix, number of dimensions, and perplexity.
	 *
	 * 	@param baseData is a matrix, observations in rows and variables in columns.
	 * 	@param outDimensions number of dimensions of tSNE space.
	 * 	@param perplexity defines the notion of neighbourhood for tSNE.
	 * 	@return TSNERunner object.
	**/
	TSNERunnerV (const arma::mat& baseData, size_t outDimensions, double perplexity);

	/**	Construct TSNERunnerV with matrix and number of dimensions
	 *
	 * 	@param baseData is a matrix, observations in rows and variables in columns.
	 * 	@param outDimensions number of dimensions of tSNE space.
	 * 	@return TSNERunner object.
	**/
	TSNERunnerV (const arma::mat& baseData, size_t outDimensions);

	~TSNERunnerV() override {};

	/**	Set number of dimensions for initial dimensionality reduction using PCA.
	 *
	 * 	@param initialDimensions number of initial dimensions. Should be less than number of rows and number of columns
	 * 		of _baseData.
	 * 	@return number of dimensions actually set (e.g. if data has lower dimensions)
	**/
	size_t setInitialDimensions (size_t initialDimensions);

	/// Get number of dimensions for initial dimensionality reduction
	size_t getInitialDimensions(void) const;

	/**	@brief Run tSNE on data
	 *
	 *	Perform tSNE dimensionality reduction on data.
	 *	@return Kullback-Leibler divergence
	**/
	double run (void) override;

	/**	Return tSNE transformation matrix.
	 *
	 * 	After the tSNE transformation has been computed with run(), this returns the matrix that transforms any
	 * 	vector (or matrix of such vectors) from the high-dimensionsional feature space to the low-dimensional
	 * 	tSNE space.
	 * 	@return An \f$n\times m\f$ matrix, where \f$n=\f$ _outDimensions and \f$m\f$ is the number of columns of _baseData.
	**/
	const arma::mat& getTSNETransform (void) const;

	/**	Project data into t-SNE space
	 *
	 * 	Project data into t-SNE space once transformation has been performed
	 * 	@param newData Matrix of data to be projected. Has to have as many
	 * 			columns as _baseData
	 * 	@return Matrix containing projected data; has as many rows as newData
	 * 			and _outDimensions columns
	**/
	arma::mat projectDataIntoTSNESpace (const arma::mat& newData) const;

protected:
	void reinitialize (void) override;

private:
	size_t _initialDimensions;
	const arma::mat& _baseData;
	arma::mat _tSNETransform; ///< Stores tSNE-transformation

};




#endif /* end of include guard: TSNERUNNERV_HPP_938F72DD */
