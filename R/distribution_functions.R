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
#' x = rvm(10, 0, 1)
#' d = dvm(x, 0, 1)
#' p = pvm(x, 0, 1)
#' @name vm
NULL

#' @rdname vm
#' @export
dvm = function(x, mu = 0, kappa = 1, log = FALSE) {
  logdens <- -log(2 * pi) - 
    log(RTMB::besselI(kappa, 0)) + 
    kappa * cos(x - mu)
  
  if(log){
    return(logdens)
  } else{
    return(exp(logdens))
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


#' wrapped Cauchy distribution
#' 
#' Density and random generation for the wrapped Cauchy distribution.
#' 
#' @details
#' The implementation of \code{dwrpcauchy} allows for automatic differentiation with \code{RTMB}. 
#' \code{rwrpcauchy} is imported from \code{CircStats}.
#'
#' @param x vector of angles measured in radians at which to evaluate the density function.
#' @param mu mean direction of the distribution measured in radians.
#' @param rho concentration parameter of the distribution, must be in the interval from 0 to 1.
#' @param log logical; if \code{TRUE}, densities are returned on the log scale.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param wrap logical; if \code{TRUE}, generated angles are wrapped to the interval [-pi, pi].
#'
#' @return \code{dwrpcauchy} gives the density and \code{rwrpcauchy} generates random deviates.
#'
#' @examples 
#' set.seed(1)
#' x = rwrpcauchy(10, 0, 1)
#' d = dwrpcauchy(x, 0, 1)
#' @name wrpcauchy
NULL

#' @rdname wrpcauchy
#' @export
dwrpcauchy <- function(x, mu = 0, rho, log = FALSE) {
  logdens <- - log(2 * pi) + 
    log(1 - rho^2) - 
    log(1 + rho^2 - 2 * rho * cos(x - mu))
  
  if(log){
    return(logdens)
  } else{
    return(exp(logdens))
  }
}

#' @rdname wrpcauchy
#' @export
#' @importFrom CircStats rwrpcauchy
rwrpcauchy = function(n, mu = 0, rho, wrap = TRUE) {
  angles = CircStats::rwrpcauchy(n, mu, rho)
  
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



#' Skew normal distribution
#' 
#' Density, distribution function, quantile function and random generation for
#' the skew normal distribution.
#' 
#' @details
#' This implementation of \code{dskewnorm} allows for automatic differentiation with \code{RTMB} while the other functions are imported from the \code{sn} package.
#'
#' @param x,q vector of quantiles
#' @param p vector of probabilities
#' @param n number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param xi location parameter
#' @param omega scale parameter, must be positive.
#' @param alpha skewness parameter, +/- \code{Inf} is allowed.
#' @param log logical; if \code{TRUE}, probabilities/ densities \eqn{p} are returned as \eqn{\log(p)}.
#' @param ... additional parameters to be passed to the \code{sn} package functions for \code{pskewnorm} and \code{qskewnorm}.
#'
#' @return
#' \code{dskewnorm} gives the density, \code{pskewnorm} gives the distribution function, \code{qskewnorm} gives the quantile function, and \code{rskewnorm} generates random deviates.
#'
#' @examples
#' x = rskewnorm(1)
#' d = dskewnorm(x)
#' p = pskewnorm(x)
#' q = qskewnorm(p)
#' @name skewnorm
NULL

#' @rdname skewnorm
#' @export
#' @importFrom RTMB dnorm
#' @importFrom RTMB pnorm
dskewnorm <- function(x, xi = 0, omega = 1, alpha = 0, log = FALSE) {
  z = (x - xi) / omega # standardised observation
  
  log_normal_density <- RTMB::dnorm(z, log = TRUE)
  log_skew_component <- log(2) - log(omega) + log(RTMB::pnorm(alpha * z))

  log_density = log_normal_density + log_skew_component
  
  if(log) {
    return(log_density)
  } else{
    return(exp(log_density))
  }
}

#' @rdname skewnorm
#' @export
#' @importFrom sn psn
pskewnorm <- function(q, xi = 0, omega = 1, alpha = 0, ...) {
  sn::psn(x = q, xi = xi, omega = omega, alpha = alpha, ...)
}

#' @rdname skewnorm
#' @export
#' @importFrom sn qsn
qskewnorm <- function(p, xi = 0, omega = 1, alpha = 0, ...) {
  sn::qsn(p = p, xi = xi, omega = omega, alpha = alpha, ...)
}

#' @rdname skewnorm
#' @export
#' @importFrom sn rsn
rskewnorm <- function(n, xi = 0, omega = 1, alpha = 0) {
  sn::rsn(n = n, xi = xi, omega = omega, alpha = alpha)
}


#' Dirichlet distribution
#' 
#' Density of the Dirichlet distribution.
#' 
#' @details
#' This implementation of \code{ddirichlet} allows for automatic differentiation with \code{RTMB}.
#'
#' @param x vector or matrix of quantiles
#' @param alpha vector or matrix of shape parameters
#' @param log logical; if \code{TRUE}, densities \eqn{p} are returned as \eqn{\log(p)}.
#'
#' @return
#' \code{ddirichlet} gives the density.
#'
#' @examples
#' ddirichlet(c(0.2, 0.3, 0.5), c(1, 2, 3))
#' @name dirichlet
NULL

#' @rdname dirichlet
#' @export
#' @import RTMB
ddirichlet <- function(x, alpha, log = TRUE) {
  # Check if x and alpha are vectors by checking if they have dimensions
  if (is.null(dim(x))) x <- matrix(x, nrow = 1)
  if (is.null(dim(alpha))) alpha <- matrix(alpha, nrow = 1)
  
  # Ensure x and alpha have the same number of columns
  if (ncol(x) != ncol(alpha)) {
    stop("x and alpha must have the same number of columns (categories).")
  }
  
  # Compute log of the multivariate beta function B(alpha) for each row
  log_B_alpha <- rowSums(lgamma(alpha)) - lgamma(rowSums(alpha))
  
  # Compute log density for each row
  log_density <- rowSums((alpha - 1) * log(x)) - log_B_alpha
  
  # Return result based on the 'log' argument
  if (log) {
    return(log_density)
  } else {
    return(exp(log_density))
  }
}



#' Zero-inflated density constructer
#' 
#' @description 
#' Constructs a zero-inflated density function from a given probability density function
#' 
#' @details
#' The definition of zero-inflation is different for discrete and continuous distributions.
#' For discrete distributions with p.m.f. \eqn{f} and zero-inflation probability \eqn{p}, we have
#' \deqn{\Pr(X = 0) = p + (1 - p) \cdot f(0),} and
#' \deqn{\Pr(X = x) = (1 - p) \cdot f(x), \quad x > 0.}
#' 
#' For continuous distributions with p.d.f. \eqn{f}, we have
#' \deqn{f_{\text{zinfl}}(x) = p \cdot \delta_0(x) + (1 - p) \cdot f(x),}
#' where \eqn{\delta_0} is the Dirac delta function at zero.
#' 
#' @param dist either a probability density function or a probability mass function
#' @param discrete logical; if \code{TRUE}, the density for \code{x = 0} will be \code{zeroprob + (1-zeroprob) * dist(0, ...)}. Otherwise it will just be \code{zeroprob}.
#' In standard cases, this will be determined automatically. For non-standard cases, set this to \code{TRUE} or \code{FALSE} depending on the type of \code{dist}. See details.
#'
#' @returns zero-inflated density function with first argument \code{x}, second argument \code{zeroprob}, and additional arguments \code{...} that will be passed to \code{dist}.
#' @export
#'
#' @examples
#' dzinorm <- zero_inflate(dnorm)
#' dzinorm(c(NA, 0, 2), 0.5, mean = 1, sd = 1)
#' 
#' zipois <- zero_inflate(dpois)
#' zipois(c(NA, 0, 1), 0.5, 1)
zero_inflate <- function(dist, discrete = NULL) {
  dist_name <- deparse(substitute(dist))
  
  # Auto-detect whether the function is discrete
  discrete_dists <- c("dpois", "dbinom", "dnbinom", "dgeom", "dmultinom", "dhyper")
  if (is.null(discrete)) {
    discrete <- dist_name %in% discrete_dists
  }
  
  # Extract function argument names, excluding `x`
  dist_args <- setdiff(names(formals(dist)), "x")
  
  # Define core function logic for both discrete and continuous cases
  core_function <- function(x, zeroprob, dist, discrete, dist_args, ...) {
    "[<-" <- RTMB::ADoverload("[<-")  # Ensure correct assignment behavior
    
    args <- list(...)
    
    # Initialize output vector
    out <- numeric(length(x))
    
    # Handle NA values properly
    is_na <- is.na(x)
    is_zero <- x == 0 & !is_na
    is_nonzero <- x != 0 & !is_na  # Exclude NA and 0
    
    # Common logic: assign probability at zero
    out[is_zero] <- zeroprob  
    
    if (discrete) {
      # Logic for discrete distributions: Handle x = 0 case with dist(0)
      out[is_zero] <- out[is_zero] + (1 - zeroprob) * do.call(dist, c(list(0), args))
    }
    
    # Logic for non-zero values
    out[is_nonzero] <- (1 - zeroprob) * do.call(dist, c(list(x[is_nonzero]), args))
    
    # Ensure NAs remain as NA
    out[is_na] <- NA
    
    return(out)
  }
  
  # Return function with core logic encapsulated
  function(x, zeroprob, ...) {
    core_function(x, zeroprob, dist, discrete, dist_args, ...)
  }
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
  
  if(is.null(dim(x_centered))){ # vector
    quadform = lambda * crossprod(x_centered, S %*% x_centered)
  } else{ # matrix
    y = S %*% t(x_centered)
    quadform = lambda * rowSums(x_centered * t(y))
  }
  
  logdens = as.numeric(0.5 * (-k * log(2*pi) + k * log(lambda) + logdetS - quadform))
  
  if(log){
    return(logdens)
  } else{
    return(exp(logdens))
  }
}
