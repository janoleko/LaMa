#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
arma::cube semigroup_cpp(const arma::mat& Q, const std::vector<double>& times) {
  int N = Q.n_cols;
  int n = times.size();
  
  arma::cube Gamma(N, N, n);
  
  for (unsigned int i = 0; i < n; i++) {
    // Assuming you want element-wise multiplication, use % for element-wise multiplication
    Gamma.slice(i) = arma::expmat(Q * times.at(i));
  }
  
  return Gamma;
}