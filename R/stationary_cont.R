#' Compute the stationary distribution of a continuous-time Markov chain
#'
#' @description
#' A well-behaved continuous-time Markov chain converges to a unique stationary distribution, here called \eqn{\pi}.
#' This distribution satisfies
#' \deqn{\pi Q = 0,} subject to \eqn{\sum_{j=1}^N \pi_j = 1},
#' where \eqn{Q} is the infinitesimal generator of the Markov chain.
#' This function solves the linear system of equations above.
#
#' @param Q infinitesimal generator matrix of dimension c(N,N)
#'
#' @return stationary distribution of the continuous-time Markov chain with given generator matrix
#' @export
#' @import RTMB
#'
#' @examples
#' Q = generator(c(-2,-2))
#' Pi = stationary_cont(Q)
stationary_cont = function(Q){
  "[<-" <- ADoverload("[<-") # overloading assignment operators, currently necessary
  "c" <- ADoverload("c")
  "diag<-" <- ADoverload("diag<-")
  
  N = dim(Q)[1]
  Pi = RTMB::solve(t(Q + 1), rep(1,N))
  names(Pi) = paste("state", 1:N)
  Pi
}