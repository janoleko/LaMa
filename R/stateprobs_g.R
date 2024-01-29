#' Calculate conditional local state probabilities for inhomogeneous HMMs
#' 
#' Computes \cr \cr
#' \eqn{\Pr(S_t = j \mid X_1, ..., X_T)} \cr \cr
#'
#' @param delta Initial distribution of length N
#' @param Gamma Array of transition probability matrices of dimension c(N,N,n-1), as in a time series of length n, there are only n-1 transitions. 
#' @param allprobs Matrix of state-dependent probabilities/ density values of dimension c(n, N)
#'
#' @return Matrix of conditional state probabilities of dimension c(n,N)
#' @export
#'
#' @examples
#' Gamma = tpm_g(runif(99), matrix(c(-1,-1,1,-2), nrow = 2, byrow = TRUE))
#' delta = c(0.5, 0.5)
#' allprobs = matrix(runif(200), nrow = 100, ncol = 2)
#' 
#' probs = stateprobs_g(delta, Gamma, allprobs)
stateprobs_g = function(delta, Gamma, allprobs){
  n = nrow(allprobs)
  N = ncol(allprobs)
  
  if(dim(Gamma)[3]==n){
    Gamma = Gamma[,,-1]
  }

  lalpha = Lcpp:::logalpha(delta, Gamma, allprobs)
  lbeta = Lcpp:::logbeta(Gamma, allprobs)
  
  c = max(lalpha[n,])
  llk = c + log(sum(exp(lalpha[n,]-c)))
  
  probs = matrix(nrow = n, ncol = N)
  for(t in 1:n){
    probs[t,] = exp(lalpha[t,] + lbeta[t,] - llk)
  }
  # rowsums should already be one
  probs = probs / rowSums(probs)
  return(probs)
}
