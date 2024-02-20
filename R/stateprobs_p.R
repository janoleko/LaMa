#' Calculate conditional local state probabilities for periodically inhomogeneous HMMs
#' 
#' Computes \cr \cr
#' \eqn{\Pr(S_t = j \mid X_1, ..., X_T)} \cr \cr
#'
#' @param delta Initial or periodically stationary distribution of length N
#' @param Gamma Array of transition probability matrices of dimension c(N,N,L) where L is the cycle length. \cr \cr
#' Here we use the definition \eqn{\Pr(S_t=j \mid S_{t-1}=i) = \gamma_{ij}^{(t)}}
#' such that the transition probabilities between time point \eqn{t-1} and \eqn{t} are an element of \eqn{\Gamma^{(t)}}.
#' @param allprobs Matrix of state-dependent probabilities/ density values of dimension c(n, N)
#' @param tod (Integer valued) time variable in 1, ..., L, mapping the data index to a generalized time of day (length n).
#' For half-hourly data L = 48. It could, however, also be day of year for daily data and L = 365.
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

stateprobs_p = function(delta, Gamma, allprobs, tod){
  n = nrow(allprobs)
  N = ncol(allprobs)
  
  Gammalarge = array(dim = c(N,N,n-1))
  
  for(t in unique(tod)){
    ind = which(tod[-1]==t)
    Gammalarge[,,ind] = Gamma[,,t]
  }
  
  lalpha = logalpha(delta, Gammalarge, allprobs)
  lbeta = logbeta(Gammalarge, allprobs)
  
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
