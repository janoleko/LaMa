// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::cube tpm_g_cpp(const arma::mat& Z, const arma::mat& beta, const arma::uword N)
{
  // get number of data points
  arma::uword n = Z.n_rows;
  
  arma::mat expEta(n, N*(N-1));
  expEta = arma::exp(Z * beta.t());
  
  // initialize matrix and cube
  arma::mat G(N,N);
  arma::cube Gamma(N,N,n);
  
  for(arma::uword i = 0; i < n; i++) {
    // Create a diagonal matrix G with N as diagonal elements
    G = arma::eye<arma::mat>(N, N);
    // Modify non-diagonal elements of G using values from expEta
    G(arma::find(G == 0)) = expEta.row(i).t();
    // Normalize each row of the matrix G by dividing it by the sum of its elements
    G = arma::normalise(G, 1, 1);
    Gamma.slice(i) = G;
  }
  return Gamma;
}
