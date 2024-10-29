#' Viterbi algorithm for state decoding in periodically inhomogeneous HMMs
#' 
#' The Viterbi algorithm allows one to decode the most probable state sequence of an HMM.
#'
#' @param delta initial distribution of length N, or matrix of dimension c(k,N) for k independent tracks, if \code{trackID} is provided
#' 
#' This could e.g. be the periodically stationary distribution (for each track).
#' @param Gamma array of transition probability matrices for each time point in the cycle of dimension c(N,N,L), where L is the length of the cycle
#' @param allprobs matrix of state-dependent probabilities/ density values of dimension c(n, N)
#' @param tod (Integer valued) variable for cycle indexing in 1, ..., L, mapping the data index to a generalised time of day (length n)
#' 
#' For half-hourly data L = 48. It could, however, also be day of year for daily data and L = 365.
#' @param trackID optional vector of k track IDs, if multiple tracks need to be decoded separately
#'
#' @return vector of decoded states of length n
#' @export
#'
#' @examples
#' delta = c(0.5, 0.5)
#' beta = matrix(c(-2, 1, -1,
#'                 -2, -1, 1), nrow = 2, byrow = TRUE)
#' Gamma = tpm_p(1:24, 24, beta)
#' 
#' tod = rep(1:24, 10)
#' n = length(tod)
#' 
#' allprobs = matrix(runif(2*n), nrow = n, ncol = 2)
#' states = viterbi_p(delta, Gamma, allprobs, tod)
viterbi_p = function(delta, Gamma, allprobs, tod, trackID = NULL){
  n = nrow(allprobs)
  N = ncol(allprobs)
  
  # creating repeating Gamma array from L unique tpms
  Gammanew = Gamma[,,tod]
  
  if(is.null(trackID)){
    Gammanew = Gammanew[,,-1]
  }
  
  # as.integer(viterbi_g_cpp(allprobs, delta, Gammanew))  
  viterbi_g(delta, Gammanew, allprobs, trackID) 
}
