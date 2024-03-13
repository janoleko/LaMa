#ifndef MY_FUNCTIONS_H
#define MY_FUNCTIONS_H

#include <Rcpp.h>
#include <armadillo>

arma::rowvec rep_times(const arma::rowvec& x, const Rcpp::IntegerVector& times);

#endif // MY_FUNCTIONS_H