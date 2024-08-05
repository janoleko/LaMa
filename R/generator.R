#' Build the generator matrix of a continuous-time Markov chain
#' 
#' @description
#' This function builds the infinitesimal generator matrix for a continuous-time Markov chain from an unconstraint parameter vector.\cr
#' 
#' Compatible with automatic differentiation by RTMB
#'
#' @param param Unconstraint parameter vector of length N*(N-1) where N is the number of states of the Markov chain
#' @param byrow Logical that indicates if the transition probability matrix should be filled by row. 
#' Defaults to FALSE, but should be set to TRUE if one wants to work with a matrix of beta parameters returned by popular HMM packages like moveHMM, momentuHMM, or hmmTMB.
#'
#' @return Infinitesimal generator matrix of dimension c(N,N)
#' @export
#' @import RTMB
#'
#' @examples
#' # 2 states: 2 free off-diagonal elements
#' generator(rep(-1, 2)
#' # 3 states: 6 free off-diagonal elements
#' generator(rep(-2, 6)
generator = function(param, byrow = FALSE) {
  
  "[<-" <- ADoverload("[<-") # currently necessary
  
  K = length(param)
  # for N > 1: N*(N-1) is bijective with solution
  N = as.integer(0.5 + sqrt(0.25 + K), 0)
  
  Q = diag(N)
  Q[!Q] = exp(param)
  diag(Q) = 0
  
  if(byrow) Q = t(Q) # transpose if necessary
  
  diag(Q) = -rowSums(Q)
  
  Q
}