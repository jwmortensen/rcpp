#include <Rcpp.h>
using namespace Rcpp;
/*
 * A few things to notice:
 *  - Vectorised functions: mean, var, etc
 *  - Don't have to worry about allocating or deallocating memory
 *  - calls to qnorm and qt (don't have to use GSL library)
 */
// [[Rcpp::export]]
int t_test_pool_rcpp(NumericVector y1, NumericVector y2, double alpha) {
  double m1 = mean(y1); 
  double m2 = mean(y2);
  double v1 = var(y1);
  double v2 = var(y2);
  int n1 = y1.size();
  int n2 = y2.size();
  int df = n1 + n2 - 2;
  double v_pool = ( (n1-1)*v1 + (n2-1)*v2 ) / df;
  double test_statistic = ( m1 - m2 ) / sqrt( v_pool * (1.0 / n1 + 1.0 / n2 ) );
  double critical_value = R::qt(1 - alpha/2, df, 1, 0); // Notice call to qt
  if (std::abs(test_statistic) >= critical_value) {
    return 1;
  } else {
    return 0;
  }
}

// [[Rcpp::export]]
int t_test_lln_rcpp(NumericVector y1, NumericVector y2, double alpha) {
  double m1 = mean(y1);
  double m2 = mean(y2);
  double v1 = var(y1);
  double v2 = var(y2);
  int n1 = y1.size();
  int n2 = y2.size();
  double test_statistic = ( m1 - m2 ) / sqrt( v1 / n1 + v2 / n2 );
  double critical_value = R::qnorm(1-alpha/2, 0, 1, 1, 0);
  if (std::abs(test_statistic) >= critical_value) {
    return 1;
  } else {
    return 0;
  }
}

// [[Rcpp::export]]
int t_test_welch_rcpp(NumericVector y1, NumericVector y2, double alpha) {
  double m1 = mean(y1);
  double m2 = mean(y2);
  double v1 = var(y1);
  double v2 = var(y2);
  int n1 = y1.size();
  int n2 = y2.size();
  double k = (v1/n1) / (v1/n1 + v2/n2);
  int df = (n1-1)*(n2-1) / ((1-k)*(1-k)*(n1-1) + k*k*(n2-1));
  double test_statistic = ( m1 - m2 ) / sqrt( v1/n1 + v2/n2 );
  double critical_value = R::qt(1 - alpha/2, df, 1, 0);
  if (std::abs(test_statistic) >= critical_value) {
    return 1;
  } else {
    return 0;
  }
}

// [[Rcpp::export]]
NumericVector engine(int n1, int n2, double mu1, double mu2, double sigma1, double sigma2, double alpha) {
  NumericVector y1(n1);
  NumericVector y2(n2);
  for (int i=0; i < n1; i++) y1[i] = mu1 + R::rnorm(0, sigma1);
  for (int i=0; i < n2; i++) y2[i] = mu2 + R::rnorm(0, sigma2);
  NumericVector result(3);
  result[0] = t_test_pool_rcpp(y1, y2, alpha);
  result[1] = t_test_lln_rcpp(y1, y2, alpha);
  result[2] = t_test_welch_rcpp(y1, y2, alpha);
  return result;
}

// [[Rcpp::export]] 
NumericMatrix simulation_rcpp(int n1, int n2, 
                         double mu1, double sigma1, 
                         NumericVector mu2, NumericVector sigma2, 
                         int nreps, double alpha) {
  int sigma2length = sigma2.size();
  int mu2length = mu2.size();
  NumericMatrix output(sigma2length*mu2length, 5);
  NumericVector result(3);
  int counter = 0;
  for (int i=0; i < sigma2length; i++) {
    for (int j=0; j < mu2length; j++) {
      double sum0 = 0;
      double sum1 = 0;
      double sum2 = 0;
      for (int k=0; k < nreps; k++) {
        result = engine(n1, n2, mu1, mu2[j], sigma1, sigma2[i], alpha);
        sum0 += result[0];
        sum1 += result[1];
        sum2 += result[2];
      }
      output(counter++, _) = NumericVector::create(mu2[j], sigma2[i], 
                                                sum0/nreps, sum1/nreps, sum2/nreps);
    }
  }
  colnames(output) = CharacterVector::create("mu2", "sigma2", "pool", "lln", "welch");
  return output;
}
