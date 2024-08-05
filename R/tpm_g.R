#' Build all transition probability matrices of an inhomogeneous HMM
#' 
#' In an HMM, we can model the influence of covariates on the state process, by linking them to the transition probabiltiy matrix. 
#' Most commonly, this is done by specifying a linear predictor \cr \cr
#' \eqn{ \eta_{ij}^{(t)} = \beta^{(ij)}_0 + \beta^{(ij)}_1 z_{t1} + \dots + \beta^{(ij)}_p z_{tp} } \cr \cr
#' for each off-diagonal element (\eqn{i \neq j}) and then applying the inverse multinomial logistic link to each row.
#' This function efficiently calculates all transition probabilty matrices for a given design matrix \eqn{Z} and parameter matrix beta.
#'
#' @param Z Covariate design matrix (excluding intercept column) of dimension c(n, p), where p can also be one (i.e. Z can be a vector).
#' @param beta Matrix of coefficients for the off-diagonal elements of the transition probability matrix.
#' Needs to be of dimension c(N*(N-1), p+1), where the first column contains the intercepts.
#' @param byrow Logical that indicates if each transition probability matrix should be filled by row. 
#' Defaults to FALSE, but should be set to TRUE if one wants to work with a matrix of beta parameters returned by popular HMM packages like moveHMM, momentuHMM, or hmmTMB.
#'
#' @return Array of transition probability matrices of dimension c(N,N,n)
#' @export
#' @import RTMB
#'
#' @examples
#' n = 1000
#' Z = matrix(runif(n*2), ncol = 2)
#' beta = matrix(c(-1, 1, 2, -2, 1, -2), nrow = 2, byrow = TRUE)
#' Gamma = tpm_g(Z, beta)
tpm_g = function(Z, beta, byrow = FALSE){
  Z = cbind(1, Z) # adding intercept column
  K = nrow(beta)
  p = ncol(beta)-1
  if(ncol(Z)!=p+1){
      stop("The dimensions of Z and beta do not match - you may have included an intercept column.")
  } else{
    # for N > 1: K = N*(N-1) is bijective with solution
    N = as.integer(0.5 + sqrt(0.25+K), 0)
    Gamma = tpm_g_cpp(Z, beta, N, byrow)
    return(Gamma)
  }
}


