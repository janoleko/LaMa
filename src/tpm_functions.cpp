// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::cube tpm_g_cpp(const arma::mat& Z, const arma::mat& beta, 
                     const arma::uword N, const bool byrow)
{
  // get number of data points
  arma::uword n = Z.n_rows;
  
  arma::mat expEta(n, N*(N-1));
  expEta = arma::exp(Z * beta.t());
  
  // initialize matrix and cube
  arma::mat G(N,N);
  arma::cube Gamma(N,N,n);
  
  if(byrow == true) {
    // if byrow is true Gamma should be created rowwise: extra transpose step
    for(arma::uword i = 0; i < n; i++) {
      // Create a diagonal matrix G with N as diagonal elements
      G = arma::eye<arma::mat>(N, N);
      // Modify non-diagonal elements of G using values from expEta
      G(arma::find(G == 0)) = expEta.row(i).t();
      // transpose for byrow
      G = G.t();
      // Normalize each row of the matrix G by dividing it by the sum of its elements
      G = arma::normalise(G, 1, 1);
      Gamma.slice(i) = G;
    }
  } else {
    for(arma::uword i = 0; i < n; i++) {
      // Create a diagonal matrix G with N as diagonal elements
      G = arma::eye<arma::mat>(N, N);
      // Modify non-diagonal elements of G using values from expEta
      G(arma::find(G == 0)) = expEta.row(i).t();
      // Normalize each row of the matrix G by dividing it by the sum of its elements
      G = arma::normalise(G, 1, 1);
      Gamma.slice(i) = G;
    }
  }
  return Gamma;
}

// [[Rcpp::export]]
arma::cube semigroup_cpp(const arma::mat& Q, const std::vector<double>& times) {
  int N = Q.n_cols;
  int n = times.size();
  
  arma::cube Gamma(N, N, n);
  
  for (unsigned int i = 0; i < n; i++) {
    Gamma.slice(i) = arma::expmat(Q * times.at(i));
  }
  
  return Gamma;
}

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
