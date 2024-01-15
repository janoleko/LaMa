#' General forward algorithm with time-varying transition probability matrix
#'
#' @param delta Initial or periodically stationary distribution (of length N)
#' @param Gamma Pre-calculated array of Gamma matrices (of dimension c(N,N,n))
#' @param allprobs allprobs matrix (of dimension c(n, N))
#'
#' @return Log-likelihood for given data and parameters
#' @export
#'
#' @examples
forward_g = function(delta, Gamma, allprobs){
  forward_cpp_g(allprobs, delta, Gamma)
}