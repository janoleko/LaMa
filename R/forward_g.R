#' R-Wrapper function for general forward algorithm with t.p.m. as an array of dim c(N,N,T)
#'
#' @param delta initial distribution
#' @param allprobs allprobs matrix
#' @param Gamma array of Gamma matrices of dim c(N,N,T)
#'
#' @return log likelihood
#' @export
#'
#' @examples
forward_g = function(delta, allprobs, Gamma){
  forward_cpp_g(allprobs, delta, Gamma)
}