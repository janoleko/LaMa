#' Build the transition probability matrix from unconstrained parameter vector
#'
#' @description
#' Markov chains are parametrised in terms of a transition probability matrix \eqn{\Gamma}, for which each row contains a conditional probability distribution of the next state given the current state.
#' Hence, each row has entries between 0 and 1 that need to sum to one. 
#' 
#' For numerical optimisation, we parametrise in terms of unconstrained parameters, thus this function computes said matrix from an unconstrained paramter vector via the inverse multinomial logistic link (also known as softmax) applied to each row.
#'
#' @param param unconstraint parameter vector of length N*(N-1) where N is the number of states of the Markov chain
#' @param byrow logical indicating if the transition probability matrix should be filled by row
#' 
#' Defaults to \code{FALSE}, but should be set to \code{TRUE} if one wants to work with a matrix of beta parameters returned by popular HMM packages like \code{moveHMM}, \code{momentuHMM}, or \code{hmmTMB}.
#'
#' @return Transition probability matrix of dimension c(N,N)
#' @export
#' @import RTMB
#'
#' @examples
#' # 2 states: 2 free off-diagonal elements
#' param1 = rep(-1, 2)
#' Gamma1 = tpm(param1)
#' 
#' # 3 states: 6 free off-diagonal elements
#' param2 = rep(-2, 6)
#' Gamma2 = tpm(param2)
tpm = function(param, byrow = FALSE) {
  
  "[<-" <- RTMB::ADoverload("[<-") # overloading assignment operators, currently necessary
  "c" <- ADoverload("c")
  "diag<-" <- ADoverload("diag<-")
  
  K = length(param)
  # for N > 1: N*(N-1) is bijective with solution
  N = as.integer(0.5 + sqrt(0.25 + K), 0)
  
  Gamma = diag(N)
  Gamma[!Gamma] = exp(param[1:(N*(N-1))])
  
  if(byrow) Gamma = t(Gamma)
  
  Gamma = Gamma / rowSums(Gamma)
  Gamma
}
