// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::sp_mat rep_times_cpp(const arma::rowvec& x, const arma::ivec& times) {
  // Validate input sizes
  if (x.n_elem == 0 || times.n_elem == 0) {
    // Handle empty vectors as needed
    return arma::sp_mat();
  }
  
  // Calculate the total size of the result vector
  arma::uword totalSize = arma::accu(arma::conv_to<arma::umat>::from(times));
  
  // Pre-allocate the result vector with the correct size
  arma::sp_mat result(1, totalSize);
  
  // Index for the result vector
  arma::uword resultIndex = 0;
  
  // Iterate through each element of 'x' and replicate according to 'times'
  for (arma::uword i = 0; i < x.n_elem; ++i) {
    for (arma::uword j = 0; j < times(i); ++j) {
      result(0, resultIndex++) = x(i);
    }
  }
  
  return result;
}