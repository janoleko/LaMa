// [[Rcpp::depends("RcppArmadillo")]]

#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
arma::rowvec rep_times(const arma::rowvec& x, const IntegerVector& times) {
  std::size_t n = times.size();
  if (n != 1 && n != x.size())
    stop("Invalid 'times' value");
  
  std::size_t n_out = std::accumulate(times.begin(), times.end(), 0);
  arma::rowvec res(n_out);
  
  auto begin = res.begin();
  for (std::size_t i = 0, ind = 0; i < n; ind += times[i], ++i) {
    auto start = begin + ind;
    auto end = start + times[i];
    std::fill(start, end, x[i]);
  }
  
  return res;
}