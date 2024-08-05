#' Build all transition probability matrices of an inhomogeneous HMM
#' 
#' In an HMM, we model the influence of covariates on the state process, by linking them to the transition probabiltiy matrix. 
#' Most commonly, this is done by specifying a linear predictor \cr \cr
#' \eqn{ \eta_{ij}^{(t)} = \beta^{(ij)}_0 + \beta^{(ij)}_1 z_{t1} + \dots + \beta^{(ij)}_p z_{tp} } \cr \cr
#' for each off-diagonal element (\eqn{i \neq j}) and then applying the inverse multinomial logistic link to each row.
#' This function efficiently calculates all transition probabilty matrices for a given design matrix \eqn{Z} and parameter matrix beta.
#'
#' @param Z Covariate design matrix with or without intercept column, i.e. of dimension c(n, p) or c(n, p+1).
#' If Z has only p columns, an intercept column of ones will be added automatically.
#' @param beta Matrix of coefficients for the off-diagonal elements of the transition probability matrix.
#' Needs to be of dimension c(N*(N-1), p+1), where the first column contains the intercepts.
#' @param byrow Logical that indicates if each transition probability matrix should be filled by row. 
#' Defaults to FALSE, but should be set to TRUE if one wants to work with a matrix of beta parameters returned by popular HMM packages like moveHMM, momentuHMM, or hmmTMB.
#' @param ad Logical, indicating whether automatic differentiation with RTMB should be used. Defaults to FALSE.
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
tpm_g = function(Z, beta, byrow = FALSE, ad = FALSE){
  
  K = nrow(beta)
  p = ncol(beta)
  N = as.integer(0.5 + sqrt(0.25 + K), 0)
  
  Z = as.matrix(Z)
  
  if(!(ncol(Z) %in% c(p, p-1))){
    stop("The dimensions of Z and beta do not match.")
  }
  
  if(ncol(Z) == p-1){
    Z = cbind(1, Z) # adding intercept column
  }
  
  if(!ad) {
    
    Gamma = tpm_g_cpp(Z, beta, N, byrow) # C++ version
    
  } else if(ad) {
    "[<-" <- ADoverload("[<-") # currently necessary
    
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

