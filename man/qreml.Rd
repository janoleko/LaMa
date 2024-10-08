% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/qreml_functions.R
\name{qreml}
\alias{qreml}
\title{Quasi restricted maximum likelihood (qREML) algorithm for models with penalized splines or simple i.i.d. random effects}
\usage{
qreml(
  pnll,
  par,
  dat,
  random,
  penalty = "lambda",
  alpha = 0,
  maxiter = 100,
  tol = 1e-05,
  control = list(reltol = 1e-10, maxit = 1000),
  silent = 1,
  saveall = FALSE
)
}
\arguments{
\item{pnll}{Penalized negative log-likelihood function that is structured as dictated by \code{RTMB} and uses the \code{\link{penalty}} function from \code{LaMa} to compute the penalty.
Needs to be a function of the named list of initial parameters \code{par} only.}

\item{par}{Named list of initial parameters. The random effects can be vectors or matrices, the latter summarising several random effects of the same structure, each one being a row in the matrix.}

\item{dat}{Initial data list that contains the data used in the likelihood function, hyperparameters, and the initial penalty strength. If the initial penalty strength vector is \strong{not} called \code{lambda}, you need to specify its name in \code{dat} using the \code{penalty} argument below.
Its length needs to match the to the total number of random effects.}

\item{random}{Vector of names of the random effects in \code{par} that are penalized.

Caution: The ordering of \code{random} needs to match the order of the random effects passed to \code{\link{penalty}} inside the likelihood function.}

\item{penalty}{Optional, name given to the penalty parameter in \code{dat}. Defaults to \code{"lambda"}.}

\item{alpha}{Optional hyperparamater for exponential smoothing of the penalty strengths. For larger values smoother convergence is to be expected but the algorithm may need more iterations. Defaults to no smoothing.}

\item{maxiter}{Maximum number of iterations in the outer optimization over the penalty strength parameters.}

\item{tol}{Convergence tolerance for the penalty strength parameters.}

\item{control}{list of control parameters for \code{\link[stats:optim]{optim}} to use in the inner optimization. Here, \code{optim} uses the BFGS method which cannot be changed.

We advise against changing the default values of \code{reltol} and \code{maxit} as this can decrease the accuracy of the Laplace approximation.}

\item{silent}{Integer silencing level: 0 corresponds to full printing of inner and outer iteratinos, 1 to printing of outer iterations only, and 2 to no printing.}

\item{saveall}{Logical, if TRUE, then all model objects from each iteration are saved in the final model object. Defaults to FALSE.}
}
\value{
Returns a model list influenced by the users report statements in \code{pnll}
}
\description{
This algorithm can be used very flexibly to fit statistical models that involve \strong{penalized splines} or simple \strong{i.i.d. random effects}, i.e. that have penalties of the form
\deqn{0.5 \sum_{i} \lambda_i b_i^T S_i b_i,}
with smoothing parameters \eqn{\lambda_i}, coefficient vectors \eqn{b_i}, and fixed penalty matrices \eqn{S_i}.\cr
The \strong{qREML} algorithm is typically much faster than REML or marginal ML using the full Laplace approximation method, but may be slightly less accurate regarding the estimation of the penalty strength parameters.
Under the hood, \code{qreml} uses the R package \code{RTMB} for automatic differentiation in the inner optimization.
The user has to specify the \strong{penalized negative log-likelihood function} \code{pnll} structured as dictated by \code{RTMB} and use the \code{\link{penalty}} function contained in \code{LaMa} to compute the quadratic-form penalty inside the likelihood.
}
\examples{
data = elephant[1:1000,] # subset

# initial parameter list
par = list(logmu = log(c(0.3, 1)), # step mean
           logsigma = log(c(0.2, 0.7)), # step sd
           beta0 = c(-2,2), # state process intercept
           betaspline = matrix(rep(0, 18), nrow = 2)) # state process spline coefs
          
# data object with initial penalty strength lambda
dat = list(step = data$step, # step length
           tod = data$tod, # time of day covariate
           N = 2, # number of states
           lambda = rep(10,2)) # initial penalty strength

# building model matrices
modmat = make_matrices(~ s(tod, bs = "cp"), 
                       data = data.frame(tod = 1:24), 
                       knots = list(tod = c(0,24))) # wrapping points
dat$Z = modmat$Z # spline design matrix
dat$S = modmat$S # penalty matrix

# penalized negative log-likelihood function
pnll = function(par) {
  getAll(par, dat) # makes everything contained available without $
  Gamma = tpm_g(Z, cbind(beta0, betaspline), ad = TRUE) # transition probabilities
  delta = stationary_p(Gamma, t = 1, ad = TRUE) # initial distribution
  mu = exp(logmu) # step mean
  sigma = exp(logsigma) # step sd
  # calculating all state-dependent densities
  allprobs = matrix(1, nrow = length(step), ncol = N)
  ind = which(!is.na(step)) # only for non-NA obs.
  for(j in 1:N) allprobs[ind,j] = dgamma2(step[ind],mu[j],sigma[j])
  -forward_g(delta, Gamma[,,tod], allprobs, ad = TRUE) +
      penalty(betaspline, S, lambda) # this does all the penalization work
}

# model fitting
mod = qreml(pnll, par, dat, random = "betaspline")
}
