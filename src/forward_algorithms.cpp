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
