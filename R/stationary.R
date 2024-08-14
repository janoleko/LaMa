#' Compute the stationary distribution of a homogeneous Markov chain
#'
#' @description
#' A homogeneous, finite state Markov chain that is irreducible and aperiodic converges to a unique stationary distribution, here called \eqn{\delta}.
#' As it is stationary, this distribution satisfies \cr \cr
#' \eqn{\delta \Gamma = \delta}, subject to \eqn{\sum_{j=1}^N \delta_j = 1}, \cr \cr
#' where \eqn{\Gamma} is the transition probability matrix. 
#' This function solves the linear system of equations above.\cr
#' 
#' Compatible with automatic differentiation by RTMB
#
#' @param Gamma Transition probability matrix of dimension c(N,N)
#'
#' @return Stationary distribution of the Markov chain with the given transition probability matrix
#' @export
#' @import RTMB
#'
#' @examples
#' Gamma = tpm(c(rep(-2,3), rep(-3,3)))
#' delta = stationary(Gamma)
stationary = function(Gamma){
  "[<-" <- RTMB::ADoverload("[<-") # overloading assignment operators, currently necessary
  "c" <- ADoverload("c")
  "diag<-" <- ADoverload("diag<-")
  
  N = dim(Gamma)[1]
  delta = RTMB::solve(t(diag(N)-Gamma+1), rep(1,N))
  names(delta) = paste("state", 1:N)
  delta
}