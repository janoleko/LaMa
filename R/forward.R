#' Forward algorithm with homogeneous transition probability matrix
#'
#' @param delta initial or stationary distribution
#' @param allprobs allprobs matrix (of dimension c(n, N))
#' @param Gamma Gamma matrix (of dimension c(N,N))
#'
#' @return log likelihood
#' @export
#'
#' @examples
forward = function(delta, allprobs, Gamma){
  forward_cpp_h(allprobs, delta, Gamma)
}