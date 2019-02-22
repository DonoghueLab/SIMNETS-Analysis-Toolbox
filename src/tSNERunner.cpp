#include "tSNEExceptions.hpp"
#include "tSNERunner.hpp"

TSNERunnerAlgoParams::TSNERunnerAlgoParams () : momSwitchIter(250), stopLyingIter(100), maxIter(700), epsilon(500.0), minGain(0.01), momentum(0.5), finalMomentum(0.8)
{
}

TSNERunner::TSNERunner ()
{
}

TSNERunner::TSNERunner (double perplexity) : _perplexity(perplexity)
{
	if (_perplexity <= 0) {
		throw badValueEx;
	}
}

TSNERunner::TSNERunner (size_t outDimensions, double perplexity) : _perplexity(perplexity), _outDimensions(outDimensions)
{
	if (_perplexity <= 0) {
		throw badValueEx;
	}
	if (_outDimensions == 0) {
		throw badValueEx;
	}
}
TSNERunner::TSNERunner (size_t outDimensions) : _outDimensions(outDimensions)
{
	if (_outDimensions == 0) {
		throw badValueEx;
	}
}

TSNERunner::TSNERunner (TSNERunnerAlgoParams params) : _params(params)
{
}

TSNERunner::TSNERunner (double perplexity, TSNERunnerAlgoParams params) : _params(params), _perplexity(perplexity)
{
	if (_perplexity <= 0) {
		throw badValueEx;
	}
}

TSNERunner::TSNERunner (size_t outDimensions, double perplexity, TSNERunnerAlgoParams params) : _params(params),  _perplexity(perplexity), _outDimensions(outDimensions)
{
	if (_perplexity <= 0) {
		throw badValueEx;
	}
	if (_outDimensions == 0) {
		throw badValueEx;
	}
}
TSNERunner::TSNERunner (size_t outDimensions, TSNERunnerAlgoParams params): _params(params), _outDimensions(outDimensions)
{
	if (_outDimensions == 0) {
		throw badValueEx;
	}
}


void TSNERunner::reinitialize(void)
{
	_ydata.reset();
	_state = Uninitialized;
	_klDiv = 0;
}

void TSNERunner::setPerplexity(double perplexity)
{
	if (perplexity > 0) {
		_perplexity = perplexity;
		reinitialize();
	}
	else {
		throw badValueEx;
	}
}

double TSNERunner::getPerplexity (void) const
{
	return _perplexity;
}

void TSNERunner::setNormalizationMethod(NormMethod preProcNorm)
{
	_preProcNorm = preProcNorm;
	reinitialize();
}

NormMethod TSNERunner::getNormalizationMethod (void) const
{
	return _preProcNorm;
}

void TSNERunner::setInitialSolution (const arma::mat& initialSolution)
{
	_ydata = initialSolution;
}


const arma::mat& TSNERunner::getTSNEData(void) const
{
	if (_state < Ready) {
		throw notCalculatedEx;
	}
	return _ydata;
}

arma::mat computeJointProbability(const arma::mat& data, NormMethod nm, double perplexity)
{
	if ((data.n_rows < 1) || (data.n_cols < 1)) {
		return arma::mat();
	}
	arma::mat X(data);
	const size_t nRows = X.n_rows;

	// normalize data
	arma::mat dd(nRows, nRows);
	if (nm == NMGlobal) {
		X -= X.min();
		X /= X.max();
	}
	else if (nm == NMByFeature) {
		X -= repmat(min(X, 0), nRows, 1);
		X /= repmat(max(X, 0) + 1e-5, nRows, 1);
	}
	for (size_t i = 0; i < nRows; ++i) {
		dd(i, i) = 0;
	}
	#pragma omp parallel for shared(dd, X)
	for (size_t i = 0; i < nRows - 1; ++i) {
		for (size_t j = i + 1; j< nRows; ++j ) {
			double dist = arma::accu(arma::pow(X.row(i) - X.row(j), 2));
			dd(i, j) = dist;
			dd(j, i) = dist;
		}
	}

	arma::mat _P = d2p(dd, perplexity, 1e-5);
	return _P;
}


arma::mat d2p(const arma::mat& D, const double u, const double tol)
{
	const size_t n = D.n_rows;
	arma::mat P(n, n, arma::fill::zeros);
	arma::vec beta(n, arma::fill::ones);
	const double logU = log(u);

	#pragma omp parallel for shared(P, beta)
	for (size_t i = 0; i < n; ++i) {
		double betamin = - std::numeric_limits<double>::infinity();
		double betamax = std::numeric_limits<double>::infinity();
		arma::uvec cols(n-1);
		if (i == 0) {
			cols = arma::linspace<arma::uvec>(i+1, n-1, n - 1) ;
		}
		else if (i == n - 1) {
			cols = arma::linspace<arma::uvec>(0, i - 1, n - 1) ;
		}
		else {
			cols = join_cols( arma::linspace<arma::uvec>(0, i-1, i) , arma::linspace<arma::uvec>(i+1, n-1, n - i - 1) );
		}
		arma::uvec rows = {i};
		std::pair<double, arma::mat> p = Hbeta(D.submat(rows, cols), beta[i]);
		double Hdiff = p.first - logU;
		size_t tries = 0;
		while ((std::abs(Hdiff) > tol) && (++tries <= 50)) {
			if (Hdiff > 0) {
				betamin = beta[i];
				if (std::isinf(betamax)) {
					beta[i] = beta[i] * 2.0;
				}
				else {
					beta[i] = (beta[i] + betamax) / 2.0;
				}
			}
			else {
				betamax = beta[i];
				if (std::isinf(betamin)) {
					beta[i] = beta[i] / 2.0;
				}
				else {
					beta[i] = (beta[i] + betamin) / 2.0;
				}
			}
			p = Hbeta(D.submat(rows, cols), beta[i]);
			Hdiff = p.first - logU;
		}
		P(rows, cols) = p.second;
	}
/*
	std::cout << "min sigma: " << min(1.0/beta)<<"\n";
	std::cout << "mean sigma: " << mean(1.0/beta)<<"\n";
	std::cout << "max sigma: " << max(1.0/beta)<<"\n";
//*/
	return P;
}


/// Function that computes the Gaussian kernel values given a vector of
/// squared Euclidean distances, and the precision of the Gaussian kernel.
/// The function also computes the perplexity of the distribution.

std::pair<double, arma::mat> Hbeta(const arma::mat& D, const double beta)
{
	arma::rowvec P(D.n_cols);
	P = exp(-D * beta);
	double sumP = fmax(arma::accu(P), 1e-8);
	double H = log(sumP) + beta * arma::accu(D % P) / sumP;
	P = P / sumP;
	return std::make_pair(H, P);

}






