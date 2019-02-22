#include "tSNEExceptions.hpp"
#include "tSNERunnerV.hpp"
#include "tSNERunnerD.hpp"

std::pair<double, arma::mat> Hbeta(const arma::mat& D, const double beta);

TSNERunnerV::TSNERunnerV (const arma::mat& baseData)
	: _baseData(baseData)
{
	setInitialDimensions(_baseData.n_rows);
}; ///< baseData is a matrix, observations in rows and variables in columns

TSNERunnerV::TSNERunnerV (const arma::mat& baseData, double perplexity)
	: TSNERunner(perplexity),
		_baseData(baseData)
{
	setInitialDimensions(_baseData.n_rows);
}

TSNERunnerV::TSNERunnerV (const arma::mat& baseData, size_t outDimensions, double perplexity)
	: TSNERunner(outDimensions, perplexity),
	 	_baseData(baseData)
{
	if (_outDimensions > _baseData.n_cols) {
		throw badValueEx;
	}
	setInitialDimensions(_baseData.n_rows);
}
TSNERunnerV::TSNERunnerV (const arma::mat& baseData, size_t outDimensions)
	: TSNERunner(outDimensions),
		_baseData(baseData)
{
	if (_outDimensions > _baseData.n_cols) {
		throw badValueEx;
	}
	setInitialDimensions(_baseData.n_rows);
}

void TSNERunnerV::reinitialize (void)
{
	_tSNETransform.reset();
	TSNERunner::reinitialize();
}

size_t TSNERunnerV::setInitialDimensions (size_t initialDimensions)
{
	reinitialize();
	return _initialDimensions = std::min<size_t>(initialDimensions, std::min<size_t>(_baseData.n_cols, _baseData.n_rows));
}

size_t TSNERunnerV::getInitialDimensions (void) const
{
	return _initialDimensions;
}


double TSNERunnerV::run(void)
{
	bool hasInitialSolution = !_ydata.empty();

	arma::mat C;
	if (_baseData.n_cols < _baseData.n_rows) {
		C = _baseData.t() * _baseData;
	}
	else {
		C = _baseData * _baseData.t() / _baseData.n_rows;
	}
	arma::vec eigval;
	arma::mat eigvec;
	arma::eig_sym(eigval, eigvec, C);

	if (any(eigval < 1e-8)) {
//		std::cout << "There are negative or very small positive eigenvalues. Results may be numerically unstable. We will try to fix this ...\n";
		eigval.elem(find(eigval < 1e-8)).fill(1e-8);
	}
	arma::uvec indices;
	try {
		indices = sort_index(eigval, "descend");
	}
	catch (...) {
		std::cout << "eigenvals:\n" << eigval << "\n";
		std::cout << "C:\n" << C << "\n";
	}

	arma::mat PCPCoeff = eigvec.cols(indices.head(_initialDimensions));
	if (!(_baseData.n_cols < _baseData.n_rows)) {
		PCPCoeff = (_baseData.t() * PCPCoeff) % arma::repmat(1.0 / sqrt(eigval(indices.head(_initialDimensions)).t() * (double)_baseData.n_rows), _baseData.n_cols, 1);
	}
	arma::mat dmatPCP = _baseData * PCPCoeff;

	TSNERunnerD tsneRunnerD(dmatPCP, _outDimensions, _perplexity, _params);
	if (hasInitialSolution) {
		tsneRunnerD.setInitialSolution(_ydata);
	}
	_klDiv = tsneRunnerD.run();
	_ydata = tsneRunnerD.getTSNEData();

	_tSNETransform = PCPCoeff * (_ydata.t() * dmatPCP * arma::pinv(dmatPCP.t() * dmatPCP)).t();

	_ydata = _baseData * _tSNETransform;

	_state = Ready;

	return _klDiv;
}


const arma::mat& TSNERunnerV::getTSNETransform (void) const
{
	if (_state < Ready) {
		throw notCalculatedEx;
	}
	return _tSNETransform;
}

arma::mat TSNERunnerV::projectDataIntoTSNESpace (const arma::mat& newData) const
{
	if (_state < Ready) {
		throw notCalculatedEx;
	}
	if (newData.n_cols != _tSNETransform.n_rows) {
		throw sizeMismatchEx;
	}
	arma::mat projection = newData * _tSNETransform;
	return projection;
}


