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

  xi = matrix(0, n, ncol = N)
  foo = delta * allprobs[1, ]
  xi[1, ] = foo / sum(foo)
  for (t in 2:n){
    foo = apply(xi[t - 1, ] * Gamma[,,tod[t]], 2, max) * allprobs[t, ]
    xi[t, ] = foo / sum(foo)
  }
  iv = numeric(n)
  iv[n] = which.max(xi[n, ])
  for (t in (n - 1):1){
    iv[t] = which.max(Gamma[, iv[t + 1], tod[t]] * xi[t, ]) 
  }
  return(iv)
}
