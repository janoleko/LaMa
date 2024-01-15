#' Forward algorithm for fruit fly application
#'
#' @param delta initial distribution
#' @param allprobs allprobs matrix (of dimension c(n, N))
#' @param Gamma1 Gamma array in LD condition (dim = c(N, N, L))
#' @param Gamma2 Gamma array in DD condition (dim = c(N, N, L))
#' @param startDD startindex of DD condition
#' @param tod time of day variable
#'
#' @return log likelihood
#' @export
#'
#' @examples
forward_flies = function(delta, allprobs, Gamma, DD, tod){
  l = forward_cpp_flies(allprobs, delta, Gamma[,,,1], Gamma[,,,2], DD, (tod-1))
  return(l)
}