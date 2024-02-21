#' Build all transition probability matrices of a periodically inhomogeneous HMM
#'
#' Given a periodically varying variable such as time of day or day of year and the associated cycle length, 
#' this function calculates the transition probability matrices by applying the inverse multinomial logistic link to linear predictors of the form \cr \cr
#' \eqn{ 
#'  \eta^{(t)}_{ij} = \beta_0^{(ij)} + \sum_{k=1}^K \bigl( \beta_{1k}^{(ij)} \sin(\frac{2 \pi k t}{L}) + \beta_{2k}^{(ij)} \cos(\frac{2 \pi k t}{L}) \bigr) } \cr \cr
#' for the off-diagonal elements (\eqn{i \neq j}).
#' This is relevant for modeling e.g. diurnal variation and the flexibility can be increased by adding smaller frequencies (i.e. increasing \eqn{K}).
#'
#' @param tod Equidistant (generalized) time of day sequence, denoting the time point in a cycle.
#' For time of day and e.g. half-hourly data, this could be 1, ..., L and L = 48, or 0.5, 1, 1.5, ..., 24 and L = 24.
#' @param L Length of one full cycle, on the scale of tod
#' @param beta Matrix of coefficients for the off-diagonal elements of the transition probability matrix.
#' Needs to be of dimension c(N*(N-1), 2*degree+1), where the first column contains the intercepts.
#' @param degree Degree of the trigonometric link function. For each additional degree, one sine and one cosine frequency are added.
#' @param Z Defaults to zero if trigonometric link should be calculated. 
#' From an efficiency perspective, this option should be used within the likelhood function, as the basis expansion should not be redundantly calculated. \cr \cr
#' Furthermore, Z can also be a pre-calculated design matrix (with p columns), when one wants to use e.g. cyclic P-splines.
#' In that case, the dimension of beta needs to be c(N*(N-1), p+1) and a penalty term should be added at the end of the negative log-likelihood.
#' @param byrow Logical that indicates if each transition probability matrix should be filled by row. 
#' Defaults to FALSE, but should be set to TRUE if one wants to work with a matrix of beta parameters returned by popular HMM packages like moveHMM, momentuHMM, or hmmTMB.
#'
#' @return Array of transition probability matrices of dimension c(N,N,length(tod))
#' @export
#'
#' @examples
#' # hourly data 
#' tod = seq(1, 24, by = 1)
#' L = 24
#' beta = matrix(c(-1, 2, -1, -2, 1, -1), nrow = 2, byrow = TRUE)
#' Gamma = tpm_p(tod, L, beta, degree = 1)
#' 
#' # half-hourly data
#' ## integer tod sequence
#' tod = seq(1, 48, by = 1)
#' L = 48
#' beta = matrix(c(-1, 2, -1, -2, 1, -1), nrow = 2, byrow = TRUE)
#' Gamma1 = tpm_p(tod, L, beta, degree = 1)
#' 
#' ## equivalent specification
#' tod = seq(0.5, 24, by = 0.5)
#' L = 24
#' beta = matrix(c(-1, 2, -1, -2, 1, -1), nrow = 2, byrow = TRUE)
#' Gamma2 = tpm_p(tod, L, beta, degree = 1)
#' 
#' Gamma1-Gamma2 # same result
#' 
#' # cubic P-splines
#' set.seed(123)
#' nk = 8 # number of basis functions
#' tod = seq(0.5, 24, by = 0.5)
#' L = 24
#' k = L * 0:nk / nk # equidistant knots
#' Z = mgcv::cSplineDes(tod, k) ## cyclic spline design matrix
#' beta = matrix(c(-1, runif(8, -2, 2), # 9 parameters per off-diagonal element
#'                  -2, runif(8, -2, 2)), nrow = 2, byrow = TRUE)
#' Gamma = tpm_p(tod, L, beta, Z = Z)
tpm_p = function(tod = 1:24, L=24, beta, degree = 1, Z = NULL, byrow = FALSE){
  K = nrow(beta)
  # for N > 1: N*(N-1) is bijective with solution
  N = as.integer(0.5 + sqrt(0.25+K), 0)
  
  if(is.null(Z)){
    Z = cbind(1, trigBasisExp(tod, L, degree))
  } else{
    Z = cbind(1, Z)
  }
  
  tpm_g_cpp(Z, beta, N, byrow)
}