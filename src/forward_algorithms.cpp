// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "my_functions.h"

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
double forward_cpp_h(const arma::mat& allprobs, const arma::rowvec& delta, const arma::mat& Gamma)
{
  int N = allprobs.n_cols;
  int nObs = allprobs.n_rows;
  arma::rowvec foo(N);
  
  foo = delta % allprobs.row(0);
  double sumfoo = sum(foo);
  double l = log(sumfoo);
  arma::rowvec phi = foo/sumfoo;
  for (unsigned int i=1; i<nObs; i++){
    foo = (phi*Gamma) % allprobs.row(i);
    sumfoo = sum(foo);
    l = l + log(sumfoo);
    phi = foo/sumfoo;
  }
  return l;
}

// [[Rcpp::export]]
double forward_cpp_g(const arma::mat& allprobs, const arma::rowvec& delta, const arma::cube& Gamma)
{
  int N = allprobs.n_cols;
  int nObs = allprobs.n_rows;
  arma::rowvec foo(N);
  
  foo = delta % allprobs.row(0);
  double sumfoo = sum(foo);
  double l = log(sumfoo);
  arma::rowvec phi = foo/sumfoo;
  for (unsigned int i=1; i<nObs; i++){
    foo = (phi*Gamma.slice(i-1)) % allprobs.row(i);
    sumfoo = sum(foo);
    l = l + log(sumfoo);
    phi = foo/sumfoo;
  }
  return l;
}

// [[Rcpp::export]]
double forward_cpp_g_tracks(const arma::mat& allprobs, const arma::mat& Delta, 
                            const arma::cube& Gamma, const IntegerVector trackInd)
{
  int N = allprobs.n_cols;
  int nObs = allprobs.n_rows;
  int K = trackInd.size();

  // forward algorithm
  double l = 0;
  arma::rowvec foo(N);
  arma::rowvec phi(N);

  unsigned int k=0; // animal index
  
  for (unsigned int t=0; t<nObs; t++){
    
    if(k<K && t==(unsigned)(trackInd(k)-1)) {
      // if 't' is the 'k'-th element of 'aInd', switch to the next track
      foo = Delta.row(k) % allprobs.row(t);
      k++;
    } else {
      foo = (phi * Gamma.slice(t)) % allprobs.row(t);
    }
    
    l = l + log(sum(foo));
    phi = foo / sum(foo);
  }
  
  return l;
}

// [[Rcpp::export]]
double forward_cpp_h_tracks(const arma::mat& allprobs, const arma::mat& Delta, 
                          const arma::cube& Gamma, const IntegerVector trackInd)
{
  int N = allprobs.n_cols;
  int nObs = allprobs.n_rows;
  int K = trackInd.size();
  
  // forward algorithm
  double l = 0;
  arma::rowvec foo(N);
  arma::rowvec phi(N);
  
  unsigned int k=0; // animal index
  
  for (unsigned int t=0; t<nObs; t++){
    
    if(k<K && t==(unsigned)(trackInd(k)-1)) {
      // if 't' is the 'k'-th element of 'aInd', switch to the next track
      foo = Delta.row(k) % allprobs.row(t);
      k++;
    } else {
      foo = (phi * Gamma.slice(k-1)) % allprobs.row(t);
    }
    
    l = l + log(sum(foo));
    phi = foo / sum(foo);
  }
  
  return l;
}

// [[Rcpp::export]]
double forward_cpp_p(const arma::mat& allprobs, const arma::rowvec& delta, const arma::cube& Gamma, const std::vector<int> tod)
{
  int N = allprobs.n_cols;
  int nObs = allprobs.n_rows;
  arma::rowvec foo(N);
  
  foo = delta % allprobs.row(0);
  double sumfoo = sum(foo);
  double l = log(sumfoo);
  arma::rowvec phi = foo/sumfoo;
  for (unsigned int i=1; i<nObs; i++){
    foo = (phi*Gamma.slice(tod[i])) % allprobs.row(i);
    sumfoo = sum(foo);
    l = l + log(sumfoo);
    
    phi = foo/sumfoo;
  }
  return l;
}

// [[Rcpp::export]]
double forward_cpp_s(const arma::mat& allprobs, const arma::rowvec& delta, const arma::mat& Gamma, const IntegerVector& agsizes)
{
  int N = allprobs.n_cols;
  int nObs = allprobs.n_rows;
  arma::rowvec foo(N);
  
  foo = delta % rep_times(allprobs.row(0), agsizes);
  double sumfoo = sum(foo);
  double l = log(sumfoo);
  arma::rowvec phi = foo/sumfoo;
  
  for (unsigned int i=1; i<nObs; i++){
    foo = (phi*Gamma) % rep_times(allprobs.row(i), agsizes);
    sumfoo = sum(foo);
    l = l + log(sumfoo);
    phi = foo/sumfoo;
  }
  return l;
}

// [[Rcpp::export]]
double forward_cpp_sp(const arma::mat& allprobs, const arma::rowvec& delta, const arma::cube& Gamma, 
                      const IntegerVector& agsizes, const std::vector<int> tod)
{
  int N = allprobs.n_cols;
  int nObs = allprobs.n_rows;
  arma::rowvec foo(N);
  
  foo = delta % rep_times(allprobs.row(0), agsizes);
  double sumfoo = sum(foo);
  double l = log(sumfoo);
  arma::rowvec phi = foo/sumfoo;
  
  for (unsigned int i=1; i<nObs; i++){
    foo = (phi*Gamma.slice(tod[i])) % rep_times(allprobs.row(i), agsizes);
    sumfoo = sum(foo);
    l = l + log(sumfoo);
    phi = foo/sumfoo;
  }
  return l;
}

// [[Rcpp::export]]
arma::mat logalpha_cpp(const arma::mat& allprobs, const arma::rowvec& delta, const arma::cube& Gamma)
{
  int N = allprobs.n_cols;
  int nObs = allprobs.n_rows;
  arma::rowvec foo(N);
  arma::mat lalpha(nObs, N);
  
  foo = delta % allprobs.row(0);
  double sumfoo = sum(foo);
  double l = log(sumfoo);
  foo = foo / sumfoo;
  lalpha.row(0) = log(foo) + l;
  
  for (unsigned int i=1; i<nObs; i++){
    foo = (foo*Gamma.slice(i-1)) % allprobs.row(i);
    sumfoo = sum(foo);
    l = l + log(sumfoo);
    foo = foo / sumfoo;
    lalpha.row(i) = log(foo) + l;
  }
  return lalpha;
}

// [[Rcpp::export]]
arma::mat logbeta_cpp(const arma::mat& allprobs, const arma::cube& Gamma) {
  int N = allprobs.n_cols;
  int nObs = allprobs.n_rows;
  arma::mat lbeta(nObs, N);
  
  arma::colvec foo = arma::ones<vec>(N) / N;
  double sumfoo = sum(foo);
  double l = log(N);
  lbeta.row(nObs - 1) = arma::zeros<rowvec>(N);
  
  for (int i = nObs - 2; i >= 0; --i) { // Loop from nObs - 1 down to 0
    foo = Gamma.slice(i) * arma::diagmat(allprobs.row(i + 1)) * foo;
    sumfoo = sum(foo);
    l = l + log(sumfoo);
    foo = foo / sumfoo;
    
    lbeta.row(i) = (l + log(foo)).t();
  }
  return lbeta;
}

// [[Rcpp::export]]
arma::colvec viterbi_g_cpp(const arma::mat& allprobs, const arma::rowvec& delta, const arma::cube& Gamma) {
  int N = allprobs.n_cols;
  int nObs = allprobs.n_rows;

  arma::mat xi(nObs, N);
  arma::rowvec foo = delta % allprobs.row(0);
  
  xi.row(0) = foo / sum(foo);
  
  for (int t = 1; t < nObs; ++t) {
    for (int j = 0; j < N; ++j) {
      foo(j) = arma::max(xi.row(t - 1).t() % Gamma.slice(t - 1).col(j)) * allprobs(t, j);
    }
    xi.row(t) = foo / arma::sum(foo);
  }
  arma::colvec iv = arma::zeros<vec>(nObs);
  iv(nObs - 1) = arma::index_max(xi.row(nObs - 1));
  
  for (int t = nObs - 2; t >= 0; --t) {
    iv(t) = arma::index_max(Gamma.slice(t).col(iv(t + 1)) % xi.row(t).t());
  }
  return(iv+1);
}