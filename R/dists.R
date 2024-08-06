#' von Mises Density Function
#' 
#' Returns the density function of the van Mises distribution evaluated at a particular value.\cr\cr
#' This implementation allows for automatic differentiation with RTMB.
#'
#' @param x vector of angles measured in radians at which to evaluate the density function.
#' @param mu mean direction of the distribution measured in radians.
#' @param kappa non-negative numeric value for the concentration parameter of the distribution.
#' @param log logical; if TRUE, densities are returned on the log scale.
#'
#' @return Returns the density function of the von Mises density distribution evaluated at x.
#' @export
#'
#' @examples 
#' x = c(0, pi/2, pi)
#' dvm(x, 0, 1)
dvm = function(x, mu, kappa, log = FALSE) {
  res = 1 / (2 * pi * RTMB::besselI(kappa, 0)) * exp(kappa * cos(x - mu))
  if(log) res = log(res)
  res
}

#' Reparametrized gamma distribution
#' 
#' Density, distribution function, quantile function and random generation for 
#' the reparametrized gamma distribution in terms of mean and standard deviation.\cr\cr
#' This implementation allows for automatic differentiation with RTMB.
#' 
#' For \code{rgamma2} defaults back to the stats implementation.
#'
#' @param x,q vector of quantiles
#' @param p vector of probabilities
#' @param n number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param mu mean parameter, must be positive scalar.
#' @param sigma standard deviation parameter, must be positive scalar.
#' @param log,log.p logical; if TRUE, probabilities/ densities \eqn{p} are returned as \eqn{\log(p)}.
#' @param lower.tail logical; if TRUE (default), probabilities are \eqn{P[X <= x]}, otherwise, \eqn{P[X > x]}.
#'
#' @return
#' \code{dgamma2} gives the density, \code{pgamma2} gives the distribution function, \code{qgamma2} gives the quantile function, and \code{rgamma2} generates random deviates.
#'
#' @examples
#' x = rgamma2(1)
#' d = dgamma2(x)
#' p = pgamma2(x)
#' q = qgamma2(p)
#' @name gamma2
NULL

#' @rdname gamma2
#' @export
dgamma2 = function(x, mu = 1, sigma = 1, log = FALSE) {
  shape = mu^2 / sigma^2
  scale = sigma^2 / mu
  RTMB::dgamma(x = x, shape = shape, scale = scale, log = log)
}

#' @rdname gamma2
#' @export
pgamma2 = function(q, mu = 1, sigma = 1, lower.tail = TRUE, log.p = FALSE) {
  shape = mu^2 / sigma^2
  scale = sigma^2 / mu
  RTMB::pgamma(q = q, shape = shape, scale = scale, lower.tail = lower.tail, log.p = log.p)
}

#' @rdname gamma2
#' @export
qgamma2 = function(p, mu = 1, sigma = 1, lower.tail = TRUE, log.p = FALSE) {
  shape = mu^2 / sigma^2
  scale = sigma^2 / mu
  RTMB::qgamma(p = p, shape = shape, scale = scale, lower.tail = lower.tail, log.p = log.p)
}

#' @rdname gamma2
#' @export
rgamma2 = function(n, mu = 1, sigma = 1) {
  shape = mu^2 / sigma^2
  scale = sigma^2 / mu
  stats::rgamma(n = n, shape = shape, scale = scale)
}