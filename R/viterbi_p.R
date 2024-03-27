#' Viterbi algorithm for decoding states of periodically inhomogeneous HMMs
#'
#' @param delta Initial or periodically statioanary distribution of length N
#' @param Gamma Array of transition probability matrices of dimension c(N,N,L), where L is the cycle length.
#' @param allprobs Matrix of state-dependent probabilities/ density values of dimension c(n, N)
#' @param tod Integer valued cyclic variable to index the transition probability matrix. 
#'
#' @return Vector of decoded states of length n
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
viterbi_p = function(delta, Gamma, allprobs, tod){
  n = nrow(allprobs)
  N = ncol(allprobs)
  
  # creating repeating Gamma array from L unique tpms
  Gammanew = array(dim = c(N,N,n-1))
  for(t in unique(tod)){
    ind = which(tod[-1]==t)
    Gammanew[,,ind] = Gamma[,,t]
  }

  as.integer(viterbi_g_cpp(allprobs, delta, Gammanew))
}
