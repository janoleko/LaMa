#' Calculate continuous time transition probabilities
#'
#' A continuous-time Markov chain is described by an infinitesimal generator matrix \eqn{Q}. 
#' When observing data at time points \eqn{t_1, \dots, t_n} the transition probabilites between \eqn{t_i} and \eqn{t_{i+1}} are caluclated as
#' \deqn{\Gamma(\Delta t_i) = \exp(Q \Delta t_i),}
#' where \eqn{\exp()} is the matrix exponential. The mapping \eqn{\Gamma(\Delta t)} is also called the \strong{Markov semigroup}.
#' This function calculates all transition matrices based on a given generator and time differences.
#'
#' @param Q infinitesimal generator matrix of the continuous-time Markov chain of dimension c(N,N)
#' @param timediff time differences between observations of length n-1 when based on n observations
#' @param ad optional logical, indicating whether automatic differentiation with \code{RTMB} should be used. By default, the function determines this itself.
#' @param report logical, indicating whether \code{Q} should be reported from the fitted model. Defaults to \code{TRUE}, but only works if \code{ad = TRUE}.
#'
#' @return array of continuous-time transition matrices of dimension c(N,N,n-1)
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
tpm_cont = function(Q, timediff, ad = NULL, report = TRUE){
  
  # report quantities for easy use later
  if(report) {
    RTMB::REPORT(Q)
  }
  
  # if ad is not explicitly provided, check if delta is an advector
  if(is.null(ad)){
    # check if Q has any of the allowed classes
    if(!any(class(Q) %in% c("advector", "numeric", "matrix", "array"))){
      stop("Q needs to be either a vector, matrix or advector.")
    }
    
    # if Q is advector, run ad version of the function
    ad = inherits(Q, "advector")
  }
  
  if(!ad) {
    
    Qube = semigroup_cpp(Q, timediff) # C++ version
    
  } else if(ad) { # ad version with RTMB
    "[<-" <- ADoverload("[<-") # currently necessary

    n = length(timediff)
    N = nrow(Q)

    Qube = array(NaN, dim = c(N, N, n))

    for(t in 1:n){
      Qube[,,t] = as.matrix(Matrix::expm(Q * timediff[t])) # Matrix::expm for AD
    }
  }
  Qube
}
