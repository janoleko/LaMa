#' Compute the periodically stationary distribution of a periodically inhomogeneous Markov chain
#'
#' @description
#' If the transition probability matrix of an inhomogeneous Markov chain varies only periodically (with period length \eqn{L}), it converges to a so-called periodically stationary distribution. 
#' This happens, because the thinned Markov chain, which has a full cycle as each time step, has homogeneous transition probability matrix
#' \deqn{\Gamma_t = \Gamma^{(t)} \Gamma^{(t+1)} \dots \Gamma^{(t+L-1)}} for all \eqn{t = 1, \dots, L.}
#' The stationary distribution for time \eqn{t} satifies \eqn{\delta^{(t)} \Gamma_t = \delta^{(t)}}.
#' 
#' This function calculates said periodically stationary distribution.
#' 
#' @param Gamma array of transition probability matrices of dimension c(N,N,L)
#' @param t integer index of the time point in the cycle, for which to calculate the stationary distribution
#' 
#' If t is not provided, the function calculates all stationary distributions for each time point in the cycle.
#' @param ad optional logical, indicating whether automatic differentiation with \code{RTMB} should be used. By default, the function determines this itself.
#'
#' @return either the periodically stationary distribution at time t or all periodically stationary distributions.
#' @export
#' @import RTMB
#'
#' @examples
# setting parameters for trigonometric link
#' L = 24
#' beta = matrix(c(-1, 2, -1, -2, 1, -1), nrow = 2, byrow = TRUE)
#' Gamma = tpm_p(1:L, L, beta, degree = 1)
#' # Periodically stationary distribution for specific time point
#' delta = stationary_p(Gamma, 4)
#'
#' # All periodically stationary distributions
#' Delta = stationary_p(Gamma)
stationary_p = function(Gamma, t = NULL, ad = NULL){
  
  N = dim(Gamma)[2]
  L = dim(Gamma)[3]
  
  # if ad is not explicitly provided, check if delta is an advector
  if(is.null(ad)){
    # check if delta has any of the allowed classes
    if(!any(class(Gamma) %in% c("advector", "array"))){
      stop("Gamma needs to be either an array or advector.")
    }
    
    # if delta is advector, run ad version of the function
    ad = inherits(Gamma, "advector")
  }
  
  if(!ad) {
    if(is.null(t)){
      Delta = matrix(nrow = L, ncol = N)
      GammaT = tpm_thinned(Gamma, 1)
      Delta[1,] = stationary(GammaT)
      for(t in 2:L){
        Delta[t,] = Delta[t-1,]%*%Gamma[,,t-1]
      }
      colnames(Delta) = paste("state", 1:N)
      return(Delta)
    } else{
      GammaT = tpm_thinned(Gamma, t)
      delta = stationary(GammaT)
      names(delta) = paste("state", 1:length(delta))
      return(delta)
    }
  } else if(ad) {
    
    "[<-" <- RTMB::ADoverload("[<-") # overloading assignment operators, currently necessary
    "c" <- ADoverload("c")
    "diag<-" <- ADoverload("diag<-")
    
    if(is.null(t)) {
      Delta = matrix(NaN, nrow = L, ncol = N)
      
      GammaT = Gamma[,,1]
      for(k in 2:L) GammaT = GammaT %*% Gamma[,,k]

      Delta[1,] = stationary(GammaT)
      
      for(t in 2:L){
        Delta[t,] = c(t(Delta[t-1,]) %*% Gamma[,,t-1])
      }
      colnames(Delta) = paste("state", 1:N)
      return(Delta)
    } else{
      GammaT = Gamma[,,t]
      if(t < L){
        for(k in (t+1):L) GammaT = GammaT %*% Gamma[,,k]
        if(t > 1){
          for(k in 1:(t-1)) GammaT = GammaT %*% Gamma[,,k]
        }
      } else if(t == L){
        for(k in 1:(L-1)) GammaT = GammaT %*% Gamma[,,k]
      }
      delta = stationary(GammaT)
      return(delta)
    }
  }
}
