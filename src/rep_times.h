#ifndef REP_TIMES_H
#define REP_TIMES_H

#include <RcppArmadillo.h>
#include <Rcpp.h>

arma::sp_mat rep_times_cpp(const arma::rowvec& x, const arma::ivec& times);

#endif  // REP_TIMES_H