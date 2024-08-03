#' von Mises Density Function
#' 
#' Returns the von Mises density function evaluated at a particular value. 
#' This implementation allows for automatic differentiation with RTMB.
#'
#' @param x vector of angles measured in radians at which to evaluate the density function.
#' @param mu mean direction of the distribution measured in radians.
#' @param kappa non-negative numeric value for the concentration parameter of the distribution.
#'
#' @return Returns the von Mises density function evaluated at x.
#' @export
#'
#' @examples 
#' x = c(0, pi/2, pi)
#' dvm(x, 0, 1)
dvm = function(x, mu, kappa) {
  1 / (2 * pi * RTMB::besselI(kappa, 0)) * exp(kappa * cos(x - mu))
}

#' Reparametrized gamma distribution
#' 
#' Returns the gamma density function reparametrized in terms of mean and standard deviation. 
#' This implementation allows for automatic differentiation with RTMB.
#'
#' @param x vector of quantiles
#' @param mu mean parameter, must be positive scalar.
#' @param sigma standard deviation parameter, must be positive scalar.
#'
#' @return Returns the gamma density function evaluated at x.
#' @export
#'
#' @examples
#' dgamma2(1, 1, 1)
dgamma2 = function(x, mu, sigma) {
  shape = mu^2 / sigma^2
  scale = sigma^2 / mu
  RTMB::dgamma(x, shape = shape, scale = scale)
}