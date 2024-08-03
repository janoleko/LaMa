#' von Mises Density Function
#' 
#' Returns the von Mises density function evaluated at a particular value. 
#' This implementation allows for automatic differentiation with RTMB.
#'
#' @param x vector of angles measured in radians at which to evaluate the density function.
#' @param mu mean direction of the distribution measured in radians.
#' @param kappa non-negative numeric value for the concentration parameter of the distribution.
#'
#' @return Returns the von Mises density function evaluated at theta.
#' @export
#'
#' @examples 
#' x = c(0, pi/2, pi)
#' dvm(x, 0, 1)
dvm = function(x, mu, kappa) {
  1 / (2 * pi * RTMB::besselI(kappa, 0)) * exp(kappa * cos(x - mu))
}