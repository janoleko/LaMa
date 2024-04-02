#' Viterbi algorithm for decoding states
#'
#' @param delta Initial or stationary distribution of length N
#' @param Gamma Transition probability matrix of dimension c(N,N)
#' @param allprobs Matrix of state-dependent probabilities/ density values of dimension c(n, N)
#'
#' @return Vector of decoded states of length n
#' @export
#'
#' @examples
#' delta = c(0.5, 0.5)
#' Gamma = matrix(c(0.9, 0.1, 0.2, 0.8), nrow = 2, byrow = TRUE)
#' allprobs = matrix(runif(200), nrow = 100, ncol = 2)
#' states = viterbi(delta, Gamma, allprobs)
viterbi = function(delta, Gamma, allprobs){
  n = nrow(allprobs)
  N = ncol(allprobs)
  Gammanew = array(Gamma, dim = c(N,N,n-1))
  as.integer(viterbi_g_cpp(allprobs, delta, Gammanew))
}



