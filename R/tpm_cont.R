#' Calculation of continuous time transition probabilities
#'
#' A continuous-time Markov chain is described by an infinitesimal generator matrix \eqn{Q}. 
#' When observing data at time points \eqn{t_1, \dots, t_n} the transition probabilites between \eqn{t_i} and \eqn{t_{i+1}} are caluclated as \cr \cr
#' \eqn{\Gamma(\Delta t_i) = \exp(Q \Delta t_i)}, \cr \cr
#' where \eqn{\exp()} is the matrix exponential. The mapping \eqn{\Gamma(\Delta t)} is also called the Markov semigroup.
#' This function calculates all transition matrices based on a given generator and time differences.
#'
#' @param Q Infinitesimal generator matrix of the continuous-time Markov chain of dimension c(N,N)
#' @param timediff Time differences between observations of length n-1 when based on n observations
#' @param ad Logical, indicating whether automatic differentiation with RTMB should be used. Defaults to FALSE.
#'
#' @return An array of transition matrices of dimension c(N,N,n-1)
#' @export
#' @import RTMB
#'
#' @examples
#' # building a Q matrix for a 3-state cont.-time Markov chain
#' Q = diag(3)
#' Q[!Q] = rexp(6)
#' diag(Q) = 0
#' diag(Q) = - rowSums(Q)
#'
#' # draw time differences
#' timediff = rexp(1000, 10)
#'
#' Gamma = tpm_cont(Q, timediff)
tpm_cont = function(Q, timediff, ad = FALSE){
  
  # if(!ad) {
    
    Qube = semigroup_cpp(Q, timediff) # C++ version
    
  # } else if(ad) { # ad version with RTMB
  #   "[<-" <- ADoverload("[<-") # currently necessary
  #   
  #   n = length(timediff)
  #   N = nrow(Q)
  #   
  #   Qube = array(NaN, dim = c(N, N, n))
  #   
  #   for(t in 1:n){
  #     Qube[,,t] = expm(Q * timediff[t]) # using RTMB::expm for AD
  #   }
  # }
  Qube
}
