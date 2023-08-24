#' R-Wrapper function for forward algorithm in CPP for periodic variation
#'
#' @param delta Initial distribution
#' @param allprobs allprobs matrix
#' @param Gamma1 Gamma matrix in LD condition (array)
#' @param Gamma2 Gamma matrix in DD condition (array)
#' @param startDD startindex of DD condition
#' @param X data
#'
#' @return
#' @export
#'
#' @examples
forward = function(delta, allprobs, Gamma1, Gamma2, startDD, X){
  l = forward_cpp(allprobs, delta, Gamma[,,,1], Gamma[,,,2], startDD, X_k$tod-1)
  return(l)
}