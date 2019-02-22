#include "tSNERunnerD.hpp"

TSNERunnerD::TSNERunnerD(const arma::mat& pw_dist)
	: _pw_dist(pw_dist)
{
}
TSNERunnerD::TSNERunnerD(const arma::mat& pw_dist, double perplexity)
	: TSNERunner(perplexity),
		_pw_dist(pw_dist)
{
}

TSNERunnerD::TSNERunnerD(const arma::mat& pw_dist, size_t outDimensions, double perplexity)
	: TSNERunner(outDimensions, perplexity),
		_pw_dist(pw_dist)
{
}

TSNERunnerD::TSNERunnerD(const arma::mat& pw_dist, size_t outDimensions)
	: TSNERunner(outDimensions),
		_pw_dist(pw_dist)
{
}

TSNERunnerD::TSNERunnerD(const arma::mat& pw_dist, TSNERunnerAlgoParams params)
	: TSNERunner(params),
		_pw_dist(pw_dist)
{
}

TSNERunnerD::TSNERunnerD(const arma::mat& pw_dist, double perplexity, TSNERunnerAlgoParams params)
	: TSNERunner(perplexity, params),
		_pw_dist(pw_dist)
{
}

TSNERunnerD::TSNERunnerD(const arma::mat& pw_dist, size_t outDimensions, double perplexity, TSNERunnerAlgoParams params)
	: TSNERunner(outDimensions, perplexity, params),
		_pw_dist(pw_dist)
{
}

TSNERunnerD::TSNERunnerD(const arma::mat& pw_dist, size_t outDimensions, TSNERunnerAlgoParams params)
	: TSNERunner(outDimensions, params),
		_pw_dist(pw_dist)
{
}


double TSNERunnerD::run (void)
{
	bool hasInitialSolution = !_ydata.empty();

	arma::mat P = computeJointProbability(_pw_dist, _preProcNorm, _perplexity);

	P += P.t();
	P *= 0.5;
	P = arma::clamp(P / ((double) P.n_rows), std::numeric_limits<double>::min(), std::numeric_limits<double>::max());

	const double klConst = arma::accu(P % log(P)); // constant in KL divergence

	if (!hasInitialSolution) {
		P *= 4;
		arma::vec eigval_P;
		arma::mat eigvec_P;
		arma::eig_sym(eigval_P, eigvec_P, P);

		if (any(eigval_P < 1e-8)) {
			//std::cout << "There are negative eigenvalues for P. Results may be numerically unstable. We will try to fix this ...\n";
			eigval_P.elem(find(eigval_P < 1e-8)).fill(1e-8);
		}
		if (eigval_P.has_nan()) {
			eigval_P.elem(find_nonfinite(eigval_P)).fill(1e-8 / 2.0);
		}
		arma::uvec indices_P;

		try {
			indices_P = sort_index(eigval_P, "descend");
		}
		catch (...) {
			std::cout << "eigenvals p:\n" << eigval_P << "\n";
			std::cout << "P:\n" << P << "\n";
			std::cout << "_pw_dist:\n" << _pw_dist << "\n";
		}


		arma::mat PCPCoeff_P = eigvec_P.cols(indices_P.head(_outDimensions));
		_ydata = P * PCPCoeff_P;
		_ydata =  _ydata / repmat(max(_ydata, 0), _ydata.n_rows, 1);

	}

	arma::mat y_incs(_ydata.n_rows, _ydata.n_cols, arma::fill::zeros);
	arma::mat gains(_ydata.n_rows, _ydata.n_cols, arma::fill::ones);
	arma::mat sum_ydata, num(_ydata.n_rows, _ydata.n_cols);
	arma::mat Q, L, y_grads;

	double currentMomentum = _params.momentum;
	for (size_t i = 0; i < _params.maxIter; ++i) {
		sum_ydata = arma::sum(arma::square(_ydata), 1);
		num = arma::repmat(sum_ydata.t(), _ydata.n_rows, 1) + arma::repmat(sum_ydata, 1, _ydata.n_rows)-2 * _ydata * _ydata.t() + 1;
		num = arma::ones<arma::mat>(num.n_rows, num.n_cols) / num;
		num.diag().fill(0.0);
		Q = arma::clamp(num / arma::accu(num), std::numeric_limits<double>::min(), std::numeric_limits<double>::max());
		L = (P - Q) % num;
		y_grads = 4* (arma::diagmat(sum(L)) - L) * _ydata;
		gains = (gains + .2) % (arma::sign(y_grads) != arma::sign(y_incs))
			+ (gains * .8) % (arma::sign(y_grads) == arma::sign(y_incs));
		gains.elem( arma::find( gains < _params.minGain)).fill( _params.minGain);
		y_incs = currentMomentum * y_incs - _params.epsilon * (gains % y_grads);
		_ydata += y_incs;
		_ydata -= arma::repmat(arma::mean(_ydata, 0), _ydata.n_rows, 1);
		if (i == _params.momSwitchIter) {
			currentMomentum = _params.finalMomentum;
		}
		if (!hasInitialSolution && (i == _params.stopLyingIter)) {
			P /= 4.0;
		}
	}

	_state = Ready;
	_klDiv = klConst - accu(P % log(Q));

	return _klDiv;
}








