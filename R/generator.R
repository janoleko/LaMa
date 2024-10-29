#' Build the generator matrix of a continuous-time Markov chain
#' 
#' This function builds the \strong{infinitesimal generator matrix} for a \strong{continuous-time Markov chain} from an unconstrained parameter vector.
#' 
#' @param param unconstrained parameter vector of length N*(N-1) where N is the number of states of the Markov chain
#' @param byrow logical indicating if the transition probability matrix should be filled by row
#' @param report logical, indicating whether the generator matrix Q should be reported from the fitted model. Defaults to \code{TRUE}, but only works if when automatic differentiation with \code{RTMB} is used.
#'
#' @return infinitesimal generator matrix of dimension c(N,N)
#' @export
#' @import RTMB
#'
#' @examples
#' # 2 states: 2 free off-diagonal elements
#' generator(rep(-1, 2))
#' # 3 states: 6 free off-diagonal elements
#' generator(rep(-2, 6))
generator = function(param, byrow = FALSE, report = TRUE) {
  
  "[<-" <- ADoverload("[<-") # currently necessary
  
  K = length(param)
  # for N > 1: N*(N-1) is bijective with solution
  N = as.integer(0.5 + sqrt(0.25 + K), 0)
  
  Q = diag(N)
  Q[!Q] = exp(param)
  diag(Q) = 0
  
  if(byrow) Q = t(Q) # transpose if necessary
  
  diag(Q) = -rowSums(Q)
  
  if(report) {
    RTMB::REPORT(Q)
  }
  
  Q
}