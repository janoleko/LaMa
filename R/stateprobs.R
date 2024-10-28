#' Calculate conditional local state probabilities for homogeneous HMMs
#' 
#' Computes \cr \cr
#' \eqn{\Pr(S_t = j \mid X_1, ..., X_T)} \cr \cr
#'
#' @param delta Initial distribution of length N, or matrix of dimension c(k,N) for k independent tracks, if \code{trackID} is provided
#' @param Gamma Transition probability matrix of dimension c(N,N) or array of transition probability matrices of dimension c(N,N,k) if \code{trackID} is provided.
#' @param allprobs Matrix of state-dependent probabilities/ density values of dimension c(n, N)
#' @param trackID Optional vector of k track IDs, if multiple tracks need to be decoded separately
#'
#' @return Matrix of conditional state probabilities of dimension c(n,N)
#' @export
#'
#' @examples
#' Gamma = tpm(c(-1,-2))
#' delta = stationary(Gamma)
#' allprobs = matrix(runif(200), nrow = 100, ncol = 2)
#' 
#' probs = stateprobs(delta, Gamma, allprobs)
stateprobs = function(delta, Gamma, allprobs, trackID = NULL){
  n = nrow(allprobs)
  N = ncol(allprobs)
  
  # inflating Gamma to use stateprobs_g
  if(is.null(trackID)){
    Gammanew = array(Gamma, dim = c(N, N, n-1))
  } else{
    uID = unique(trackID)
    k = length(uID) # number of tracks
    
    if(dim(Gamma)[3] != k) stop("Number of distinct transition matrices does not match the number of tracks.")
    
    ## construct integer trackID
    integer_ID = match(trackID, uID)
    
    Gammanew = Gamma[,,integerID]
  }
  
  stateprobs_g(delta, Gammanew, allprobs, trackID)
}