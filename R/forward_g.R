#' General forward algorithm with time-varying transition probability matrix
#'
#' @param delta initial or periodically stationary distribution
#' @param allprobs allprobs matrix (of dimension c(n, N))
#' @param Gamma array of Gamma matrices (of dimension c(N,N,n))
#'
#' @return log likelihood
#' @export
#'
#' @examples
forward_g = function(delta, allprobs, Gamma){
  forward_cpp_g(allprobs, delta, Gamma)
}