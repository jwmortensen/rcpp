// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

arma::mat mvrnorm_arma(int n, arma::vec mu, arma::mat sigma) {
  int ncols = sigma.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);
}

// [[Rcpp::export]]
List gibbs_sampler_rcpp(arma::mat X, arma::vec y, int ndraws, double beta_s2) {
  int n_beta = X.n_cols;
  int n = y.size();
  arma::mat beta_draws(ndraws, n_beta, arma::fill::zeros);
  arma::vec sigma2_draws(ndraws, arma::fill::ones);
  arma::mat Lambda(n_beta, n_beta);
  arma::mat Sigma(n_beta, n_beta);
  arma::vec mu(n_beta);
  for (int i=1; i < ndraws; i++) {
    double s2_inv = 1 / sigma2_draws[i-1];
    Lambda = s2_inv * X.t() * X + (1.0 / beta_s2);
    Sigma = arma::inv(Lambda);
    mu = s2_inv * Sigma * X.t() * y;
    arma::vec beta = arma::vec(mvrnorm_arma(1, mu, Sigma).t());
    beta_draws.row(i) = beta.t();
    double shape = 2.001 + (double)  n / 2.0;
    double rate = arma::as_scalar(1.001 + 0.5 * (y - X * beta).t() * (y - X * beta));
    sigma2_draws(i) = 1.0 / R::rgamma(shape, 1.0 / rate);
  }

  return List::create(_["beta_draws"] = beta_draws, _["sigma2_draws"] = sigma2_draws);
}
