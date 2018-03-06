// [[Rcpp::depends(RcppGSL)]]
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <RcppGSL.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
using namespace Rcpp;


// [[Rcpp::export]]
arma::mat mvrnorm_arma(int n, arma::vec mu, arma::mat sigma) {
  int ncols = sigma.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);
}


RcppGSL::Vector single_mvrnorm_gsl(gsl_rng *r, RcppGSL::Vector mu, const gsl_matrix *chol) {
  int n = mu.size();
  RcppGSL::Vector Y(n);
  for (int i=0; i < n; i++) {
    Y[i] = gsl_ran_ugaussian(r);
  }
  gsl_blas_dtrmv(CblasLower, CblasNoTrans, CblasNonUnit, chol, Y);
  gsl_vector_add(Y, mu);
  return Y;
}

// [[Rcpp::export]]
RcppGSL::Matrix mvrnorm_gsl(int n, RcppGSL::Vector mu, RcppGSL::Matrix sigma) {
  int ncol = sigma.ncol();
  gsl_rng * r = gsl_rng_alloc(gsl_rng_taus);
  gsl_matrix *chol = gsl_matrix_alloc(ncol, ncol);
  gsl_matrix_memcpy(chol, sigma);
  gsl_linalg_cholesky_decomp(chol);
  RcppGSL::Matrix out(n, ncol);
  for (int i=0; i < n; i++) {
    RcppGSL::Vector vec = single_mvrnorm_gsl(r, mu, chol);
    for (int j = 0; j < vec.size(); j++) {
      out(i, j) = vec[j];
    }
  }
  free(chol);
  return out;
}