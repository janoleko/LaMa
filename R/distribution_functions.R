#' von Mises distribution
#' 
#' Density, distribution function and random generation for the von Mises distribution.
#' 
#' @details
#' The implementation of \code{dvm} allows for automatic differentiation with \code{RTMB}. 
#' \code{rvm} and \code{pvm} are imported from \code{CircStats} and \code{circular} respectively.
#'
#' @param x,q vector of angles measured in radians at which to evaluate the density function.
#' @param mu mean direction of the distribution measured in radians.
#' @param kappa non-negative numeric value for the concentration parameter of the distribution.
#' @param log logical; if \code{TRUE}, densities are returned on the log scale.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param tol the precision in evaluating the distribution function
#' @param from value from which the integration for CDF starts. If \code{NULL}, is set to \code{mu - pi}.
#' @param wrap logical; if \code{TRUE}, generated angles are wrapped to the interval [-pi, pi].
#'
#' @return \code{dvm} gives the density, \code{pvm} gives the distribution function, and \code{rvm} generates random deviates.
#'
#' @examples 
#' set.seed(1)
#' x = rvm(1000, 0, 1)
#' d = dvm(x, 0, 1)
#' p = pvm(x, 0, 1)
#' @name vm
NULL

#' @rdname vm
#' @export
dvm = function(x, mu = 0, kappa = 1, log = FALSE) {
  res = 1 / (2 * pi * RTMB::besselI(kappa, 0)) * exp(kappa * cos(x - mu))
  if(log){
    return(log(res))
  } else{
    return(res)
  }
}

#' @rdname vm
#' @export
#' @importFrom circular pvonmises
pvm = function(q, mu = 0, kappa = 1, from = NULL, tol = 1e-20) {
  # NA handling
  ind = which(!is.na(q))
  
  if(is.matrix(mu)){
    mu = mu[ind,]
  }
  if(is.matrix(kappa)){
    kappa = kappa[ind,]
  }
  
  probs = numeric(length(q))
  
  suppressWarnings(
    probs[ind] <- circular::pvonmises(q[ind], mu, kappa, from = from, tol = tol)
  )
  
  probs[-ind] = NA
  
  as.numeric(probs)
}

#' @rdname vm
#' @export
#' @importFrom CircStats rvm
rvm = function(n, mu = 0, kappa = 1, wrap = TRUE) {
  angles = CircStats::rvm(n, mu, kappa)
  
  # if generated angels should be wrapped, i.e. mapped to interval [-pi, pi], do so
  if(wrap){
    angles = (angles + pi) %% (2 * pi) - pi
  }
  angles
}


#' Reparametrised gamma distribution
#' 
#' Density, distribution function, quantile function and random generation for
#' the gamma distribution reparametrised in terms of mean and standard deviation.
#' 
#' @details
#' This implementation allows for automatic differentiation with \code{RTMB}.
#'
#' @param x,q vector of quantiles
#' @param p vector of probabilities
#' @param n number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param mean mean parameter, must be positive scalar.
#' @param sd standard deviation parameter, must be positive scalar.
#' @param log,log.p logical; if \code{TRUE}, probabilities/ densities \eqn{p} are returned as \eqn{\log(p)}.
#' @param lower.tail logical; if \code{TRUE}, probabilities are \eqn{P[X <= x]}, otherwise, \eqn{P[X > x]}.
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
#' @importFrom RTMB dgamma
dgamma2 = function(x, mean = 1, sd = 1, log = FALSE) {
  shape = mean^2 / sd^2
  scale = sd^2 / mean
  RTMB::dgamma(x = x, shape = shape, scale = scale, log = log)
}

#' @rdname gamma2
#' @export
pgamma2 = function(q, mean = 1, sd = 1, lower.tail = TRUE, log.p = FALSE) {
  shape = mean^2 / sd^2
  scale = sd^2 / mean
  RTMB::pgamma(q = q, shape = shape, scale = scale, lower.tail = lower.tail, log.p = log.p)
}

#' @rdname gamma2
#' @export
qgamma2 = function(p, mean = 1, sd = 1, lower.tail = TRUE, log.p = FALSE) {
  shape = mean^2 / sd^2
  scale = sd^2 / mean
  RTMB::qgamma(p = p, shape = shape, scale = scale, lower.tail = lower.tail, log.p = log.p)
}

#' @rdname gamma2
#' @export
rgamma2 = function(n, mean = 1, sd = 1) {
  shape = mean^2 / sd^2
  scale = sd^2 / mean
  stats::rgamma(n = n, shape = shape, scale = scale)
}


#' Reparametrised multivariate Gaussian distribution
#'
#' Density function of the multivariate Gaussian distribution reparametrised in terms of its precision matrix (inverse variance).
#' This implementation is particularly useful for defining the \strong{joint log-likelihood} with penalised splines or i.i.d. random effects that have a multivariate Gaussian distribution with fixed precision/ penalty matrix \eqn{\lambda S}.
#' As \eqn{S} is fixed and only scaled by \eqn{\lambda}, it is more efficient to precompute the determinant of \eqn{S} (for the normalisation constant) and only scale the quadratic form by \eqn{\lambda}
#' when multiple spline parameters/ random effects with different \eqn{\lambda}'s but the same penalty matrix \eqn{S} are evaluated.
#'
#' This implementation allows for automatic differentiation with \code{RTMB}.
#'
#' @param x density evaluation point, either a vector or a matrix
#' @param mu mean parameter. Either scalar or vector
#' @param S unscaled precision matrix
#' @param lambda precision scaling parameter
#' 
#' Can be a vector if \code{x} is a matrix. Then each row of \code{x} is evaluated with the corresponding \code{lambda}.
#' This is benefitial from an efficiency perspective because the determinant of \code{S} is only computed once.
#' @param logdetS Optional precomputed log determinant of the precision matrix \code{S}. If the precision matrix does not depend on parameters, it can be precomputed and passed to the function.
#' @param log logical; if \code{TRUE}, densities are returned on the log scale.
#'
#' @return vector of density values
#' @export
#' @import RTMB
#'
#' @examples
#' x = matrix(runif(30), nrow = 3)
#'
#' # iid random effects
#' S = diag(10)
#' sigma = c(1, 2, 3) # random effect standard deviations
#' lambda = 1 / sigma^2
#' d = dgmrf2(x, 0, S, lambda)
#'
#' # P-splines
#' L = diff(diag(10), diff = 2) # second-order difference matrix
#' S = t(L) %*% L
#' lambda = c(1,2,3)
#' d = dgmrf2(x, 0, S, lambda, log = TRUE)
dgmrf2 = function(x, 
                  mu = 0, 
                  S, 
                  lambda, 
                  logdetS = NULL,
                  log = FALSE) {
  # if logdet not specified, compute generalised determinant
  if(is.null(logdetS)){
    logdetS = gdeterminant(S)
  }
  k = nrow(S)
  
  x_centered = x - mu # center data
  
  if(is.matrix(x_centered)){
    y = S %*% t(x_centered)
    quadform = lambda * rowSums(x_centered * t(y))
  } else{
    quadform = lambda * crossprod(x_centered, S %*% x_centered)
  }
  
  logdens = as.numeric(0.5 * (-k * log(2*pi) + k * log(lambda) + logdetS - quadform))
  
  if(log){
    return(logdens)
  } else{
    return(exp(logdens))
  }
}
