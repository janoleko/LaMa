#' Calculation of continuous time transition probabilities
#'
#' A continuous-time Markov chain is described by an infinitesimal generator matrix \eqn{Q}. 
#' When observing data at time points \eqn{t_1, \dots, t_n} the transition probabilites between \eqn{t_i} and \eqn{t_{i+1}} are caluclated as \cr \cr
#' \eqn{\Gamma(\Delta t_i) = \exp(Q \Delta t_i)}, \cr \cr
#' where \eqn{\exp()} is the matrix exponential. The mapping \eqn{\Gamma(\Delta t)} is also called the Markov semigroup.
#' This function calculated all transition matrices based on a given generator and time differences.
#'
#' @param Q Infinitesimal generator matrix of the continuous-time Markov chain of dimension c(N,N)
#' @param times Time differences between observations of length n-1 when based on n observations
#'
#' @return An array of transition matrices of dimension c(N,N,n-1)
#' @export
#'
#' @examples
#' # building a Q matrix for a 3-state cont.-time Markov chain
#' Q = diag(3)
#' Q[!Q] = rexp(6)
#' diag(Q) = 0
#' diag(Q) = - rowSums(Q)
#'
#' # draw time differences
#' times = rexp(1000, 10)
#'
#' Gamma = semigroup(Q, times)
semigroup = function(Q, times){
  semigroup_cpp(Q, times)
}
