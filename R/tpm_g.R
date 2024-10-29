#' Build all transition probability matrices of an inhomogeneous HMM
#' 
#' In an HMM, we model the influence of covariates on the state process by linking them to the transition probabiltiy matrix. 
#' Most commonly, this is done by specifying a linear predictor
#' \deqn{ \eta_{ij}^{(t)} = \beta^{(ij)}_0 + \beta^{(ij)}_1 z_{t1} + \dots + \beta^{(ij)}_p z_{tp} }
#' for each off-diagonal element (\eqn{i \neq j}) of the transition probability matrix and then applying the inverse multinomial logistic link (also known as softmax) to each row.
#' This function efficiently calculates all transition probabilty matrices for a given design matrix \eqn{Z} and parameter matrix beta.
#'
#' @param Z covariate design matrix with or without intercept column, i.e. of dimension c(n, p) or c(n, p+1)
#' 
#' If Z has only p columns, an intercept column of ones will be added automatically.
#' @param beta matrix of coefficients for the off-diagonal elements of the transition probability matrix
#' 
#' Needs to be of dimension c(N*(N-1), p+1), where the first column contains the intercepts.
#' @param byrow logical indicating if each transition probability matrix should be filled by row
#'  
#' Defaults to \code{FALSE}, but should be set to \code{TRUE} if one wants to work with a matrix of beta parameters returned by popular HMM packages like \code{moveHMM}, \code{momentuHMM}, or \code{hmmTMB}.
#' @param ad optional logical, indicating whether automatic differentiation with \code{RTMB} should be used. By default, the function determines this itself.
#' @param report logical, indicating whether the coefficient matrix \code{beta} should be reported from the fitted model. Defaults to \code{TRUE}, but only works if \code{ad = TRUE}.
#'
#' @return array of transition probability matrices of dimension c(N,N,n)
#' @export
#' @import RTMB
#'
#' @examples
#' n = 1000
#' Z = matrix(runif(n*2), ncol = 2)
#' beta = matrix(c(-1, 1, 2, -2, 1, -2), nrow = 2, byrow = TRUE)
#' Gamma = tpm_g(Z, beta)
tpm_g = function(Z, beta, byrow = FALSE, ad = NULL, report = TRUE){
  
  K = nrow(beta)
  p = ncol(beta) - 1
  N = as.integer(0.5 + sqrt(0.25 + K), 0)
  
  Z = as.matrix(Z)
  
  if(ncol(Z) == p){
    Z = cbind(1, Z) # adding intercept column
  } else if(ncol(Z) != p + 1){
    stop("The dimensions of Z and beta do not match.")
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
  
  if(!ad) {

    Gamma = tpm_g_cpp(Z, beta, N, byrow) # C++ version
    
  } else if(ad) {
    "[<-" <- ADoverload("[<-") # overloading assignment operators, currently necessary
    "c" <- ADoverload("c")
    "diag<-" <- ADoverload("diag<-")
    
    # if(report) {
    #   RTMB::REPORT(beta) # reporting coefficient matrix
    # }
    
    expEta = exp(Z %*% t(beta))
    Gamma = array(NaN, dim = c(N, N, nrow(expEta)))
    
    for(t in 1:nrow(expEta)){
      G = diag(N)
      G[!G] = expEta[t,]
      
      if(byrow) G = t(G) # transpose if necessary
      
      Gamma[,,t] = G / rowSums(G)
    }
  }
  
  Gamma
}

