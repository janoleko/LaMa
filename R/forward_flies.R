#' Forward algorithm for fruit fly application
#'
#' @param delta Initial or periodically stationary distribution of length N
#' @param Gamma Array of transition probability matrices of dimension c(N,N,L,2).
#' @param allprobs Matrix of state-dependent probabilities/ density values of dimension c(n, N)
#' @param DD Index where DD condition starts
#' @param tod time of day variable
#'
#' @return Log-likelihood for given data and parameters
#' @export
#'
#' @examples
forward_flies = function(delta, allprobs, Gamma, DD, tod){
  l = forward_cpp_flies(allprobs, delta, Gamma[,,,1], Gamma[,,,2], DD, (tod-1))
  return(l)
}