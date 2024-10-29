#' Build all transition probability matrices of a periodically inhomogeneous HMM
#'
#' @description
#' Given a periodically varying variable such as time of day or day of year and the associated cycle length, 
#' this function calculates the transition probability matrices by applying the inverse multinomial logistic link (also known as softmax) to linear predictors of the form
#' \deqn{ 
#'  \eta^{(t)}_{ij} = \beta_0^{(ij)} + \sum_{k=1}^K \bigl( \beta_{1k}^{(ij)} \sin(\frac{2 \pi k t}{L}) + \beta_{2k}^{(ij)} \cos(\frac{2 \pi k t}{L}) \bigr) }
#' for the off-diagonal elements (\eqn{i \neq j}) of the transition probability matrix.
#' This is relevant for modeling e.g. diurnal variation and the flexibility can be increased by adding smaller frequencies (i.e. increasing \eqn{K}).
#' 
#' @details
#' Note that using this function inside the negative log-likelihood function is convenient, but it performs the basis expansion into sine and cosine terms each time it is called. 
#' As these do not change during the optimisation, using \code{\link{tpm_g}} with a pre-calculated (by \code{\link{trigBasisExp}}) design matrix would be more efficient.
#'
#' @param tod equidistant sequence of a cyclic variable
#' 
#' For time of day and e.g. half-hourly data, this could be 1, ..., L and L = 48, or 0.5, 1, 1.5, ..., 24 and L = 24.
#' @param L length of one full cycle, on the scale of tod
#' 
#' @param beta matrix of coefficients for the off-diagonal elements of the transition probability matrix
#' 
#' Needs to be of dimension c(N *(N-1), 2*degree+1), where the first column contains the intercepts.
#' @param degree degree of the trigonometric link function
#' 
#' For each additional degree, one sine and one cosine frequency are added.
#' @param Z pre-calculated design matrix (excluding intercept column)
#' 
#' Defaults to \code{NULL} if trigonometric link should be calculated. 
#' From an efficiency perspective, Z should be pre-calculated within the likelihood function, as the basis expansion should not be redundantly calculated. This can be done by using \code{\link{trigBasisExp}}.
#' @param byrow logical indicating if each transition probability matrix should be filled by row
#'  
#' Defaults to \code{FALSE}, but should be set to \code{TRUE} if one wants to work with a matrix of \code{beta} parameters returned by popular HMM packages like \code{moveHMM}, \code{momentuHMM}, or \code{hmmTMB}.
#' @param ad optional logical, indicating whether automatic differentiation with RTMB should be used. By default, the function determines this itself.
#' @param report logical, indicating whether the coefficient matrix \code{beta} should be reported from the fitted model. Defaults to \code{TRUE}, but only works if \code{ad = TRUE}.
#'
#' @return array of transition probability matrices of dimension c(N,N,length(tod))
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
tpm_p = function(tod = 1:24, L=24, beta, degree = 1, Z = NULL, byrow = FALSE, ad = NULL, report = TRUE){
  K = nrow(beta)
  p = ncol(beta) - 1 # number of covariates
  # for N > 1: K = N*(N-1) is bijective with solution
  N = as.integer(0.5 + sqrt(0.25+K), 0)
  
  # check wheter design matrix is provided
  if(is.null(Z)){ # if not, build trigonometric design matrix
    Z = cbind(1, trigBasisExp(tod, L, degree))
  } else{
    if(ncol(Z) == p){ # intercept column is missing
      Z = cbind(1, Z) # adding intercept column
    } else if(ncol(Z) != p + 1){
      stop("The dimensions of Z and beta do not match.")
    }
  }
  
  # report quantities for easy use later
  if(report) {
    RTMB::REPORT(beta)
  }
  
  # if ad is not explicitly provided, check if delta is an advector
  if(is.null(ad)){
    # check if delta has any of the allowed classes
    if(!any(class(beta) %in% c("advector", "numeric", "matrix", "array"))){
      stop("beta needs to be either a matrix or advector.")
    }
    
    # if delta is advector, run ad version of the function
    ad = inherits(beta, "advector")
  }
  
  if(!ad){
    Gamma = tpm_g_cpp(Z, beta, N, byrow)
  } else if(ad) {
    "[<-" <- ADoverload("[<-") # overloading assignment operators, currently necessary
    "c" <- ADoverload("c")
    "diag<-" <- ADoverload("diag<-")
    
    Gamma = tpm_g(Z, beta, byrow, ad = TRUE, report = FALSE)
  }
  
  Gamma
}