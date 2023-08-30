#' R-Wrapper function for forward algorithm with homogeneous t.p.m.
#'
#' @param delta initial distribution
#' @param allprobs allprobs matrix
#' @param Gamma1 Gamma matrix
#'
#' @return log likelihood
#' @export
#'
#' @examples
forward_h = function(delta, allprobs, Gamma){
  forward_cpp_h(allprobs, delta, Gamma)
}