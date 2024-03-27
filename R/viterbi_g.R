#' Viterbi algorithm for decoding states of inhomogeneous HMMs
#'
#' @param delta Initial distribution of length N
#' @param Gamma Array of transition probability matrices of dimension c(N,N,n-1), as in a time series of length n, there are only n-1 transitions. 
#' If you provide an array of dimension c(N,N,n), the first slice will be ignored. \cr
#' @param allprobs Matrix of state-dependent probabilities/ density values of dimension c(n, N)
#'
#' @return Vector of decoded states of length n
#' @export
#'
#' @examples
#' delta = c(0.5, 0.5)
#' Gamma = array(dim = c(2,2,99))
#' for(t in 1:99){
#'   gammas = rbeta(2, shape1 = 0.4, shape2 = 1)
#'   Gamma[,,t] = matrix(c(1-gammas[1], gammas[1], 
#'                       gammas[2], 1-gammas[2]), nrow = 2, byrow = TRUE)
#' }
#' allprobs = matrix(runif(200), nrow = 100, ncol = 2)
#' states = viterbi_g(delta, Gamma, allprobs)
viterbi_g = function(delta, Gamma, allprobs){
  n = nrow(allprobs)
  N = ncol(allprobs)
  
  if(dim(Gamma)[3]==n){
    # not using the first slice of Gamma, if n slices are provided
    Gamma = Gamma[,,-1]
  }
  
  as.integer(viterbi_g_cpp(allprobs, delta, Gamma))
}
