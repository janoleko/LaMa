#' Calculate conditional local state probabilities for periodically inhomogeneous HMMs
#' 
#' Computes \cr \cr
#' \eqn{\Pr(S_t = j \mid X_1, ..., X_T)} \cr \cr
#'
#' @param delta Initial distribution of length N, or matrix of dimension c(k,N) for k independent tracks, if \code{trackID} is provided. This could e.g. be the periodically stationary distribution (for each track).
#' @param Gamma Array of transition probability matrices for each time point in the cycle of dimension c(N,N,L), where L is the length of the cycle.
#' @param allprobs Matrix of state-dependent probabilities/ density values of dimension c(n, N)
#' @param tod (Integer valued) time variable in 1, ..., L, mapping the data index to a generalized time of day (length n).
#' For half-hourly data L = 48. It could, however, also be day of year for daily data and L = 365.
#' @param trackID Optional vector of k track IDs, if multiple tracks need to be decoded separately
#' 
#' @return Matrix of conditional state probabilities of dimension c(n,N)
#' @export
#'
#' @examples
#' L = 24
#' beta = matrix(c(-1, 2, -1, -2, 1, -1), nrow = 2, byrow = TRUE)
#' Gamma = tpm_p(1:L, L, beta, degree = 1)
#' delta = stationary_p(Gamma, 1)
#' allprobs = matrix(runif(200), nrow = 100, ncol = 2)
#' tod = rep(1:24, 5)[1:100]
#' 
#' probs = stateprobs_p(delta, Gamma, allprobs, tod)

stateprobs_p = function(delta, Gamma, allprobs, tod, trackID = NULL){
  n = nrow(allprobs)
  N = ncol(allprobs)
  
  Gammanew = Gamma[,,tod] # select the transition matrix for the current time of day
  
  if(is.null(trackID)){
    Gammanew = Gammanew[,,-1]
  }
  
  stateprobs_g(delta, Gammanew, allprobs, trackID)
}
