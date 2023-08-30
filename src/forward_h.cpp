#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double forward_cpp_h(arma::mat allprobs, arma::rowvec delta, arma::mat Gamma)
{
  
  int N = allprobs.n_cols;
  int nObs = allprobs.n_rows;
  arma::rowvec foo(N);
  
  foo = delta % allprobs.row(0);
  double mllk_scale = log(sum(foo));
  arma::rowvec phi = foo/sum(foo);
  for (unsigned int i=1; i<nObs; i++){
    foo = (phi*Gamma) % allprobs.row(i);
    mllk_scale = mllk_scale + log(sum(foo));
    phi = foo/sum(foo);
  }
  return mllk_scale;
}
