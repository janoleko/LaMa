// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::cube semigroup_sym_cpp(const arma::mat& Q, const std::vector<double>& times) {
  int N = Q.n_cols;
  int n = times.size();
  
  arma::cube Gamma(N, N, n);
  
  for (unsigned int i = 0; i < n; i++) {
    Gamma.slice(i) = arma::expmat_sym(Q * times.at(i));
  }
  
  return Gamma;
}