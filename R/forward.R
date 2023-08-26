#' R-Wrapper function for forward algorithm in CPP for periodic variation
#'
#' @param delta Initial distribution
#' @param allprobs allprobs matrix
#' @param Gamma1 Gamma array in LD condition (dim = c(N, N, L))
#' @param Gamma2 Gamma array in DD condition (dim = c(N, N, L))
#' @param startDD startindex of DD condition
#' @param tod time of day variable
#'
#' @return log likelihood
#' @export
#'
#' @examples
forward = function(delta, allprobs, Gamma1, Gamma2, DD, tod){
  l = forward_cpp(allprobs, delta, Gamma[,,,1], Gamma[,,,2], DD, (tod-1))
  return(l)
}