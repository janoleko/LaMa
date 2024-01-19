// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::mat tpm_thinned_t_cpp(const arma::cube& Gamma, const int t)
{
  arma::uword L = Gamma.n_slices;
  arma::uword N = Gamma.n_rows;
  
  arma::mat GammaT(N,N); 
  // first half of the product form t to L
  GammaT = Gamma.slice(t-1);
  for(arma::uword i=t; i<L; i++) {
    GammaT = GammaT * Gamma.slice(i);
  }
  // second half of the product from 1 to t-1
  for(arma::uword i=0; i<t-1; i++) {
    GammaT = GammaT * Gamma.slice(i);
  }
  return GammaT;
}
