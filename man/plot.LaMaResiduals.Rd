% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/residual_functions.R
\name{plot.LaMaResiduals}
\alias{plot.LaMaResiduals}
\title{Plot pseudo-residuals}
\usage{
\method{plot}{LaMaResiduals}(x, hist = FALSE, col = "darkblue", lwd = 1.5, main = NULL, ...)
}
\arguments{
\item{x}{pseudo-residuals as returned by \code{\link{pseudo_res}}}

\item{hist}{logical, if \code{TRUE}, adds a histogram of the pseudo-residuals}

\item{col}{character, color for the QQ-line (and density curve if \code{histogram = TRUE})}

\item{lwd}{numeric, line width for the QQ-line (and density curve if \code{histogram = TRUE})}

\item{main}{optional character vector of main titles for the plots of length 2 (or 3 if \code{histogram = TRUE})}

\item{...}{currently ignored. For method consistency}
}
\value{
NULL, plots the pseudo-residuals in a 2- or 3-panel layout
}
\description{
Plot pseudo-residuals computed by \code{\link{pseudo_res}}.
}
\examples{
## pseudo-residuals for the trex data
step = trex$step[1:200]

nll = function(par){
  getAll(par)
  Gamma = tpm(logitGamma)
  delta = stationary(Gamma)
  mu = exp(logMu); REPORT(mu)
  sigma = exp(logSigma); REPORT(sigma)
  allprobs = matrix(1, length(step), 2)
  ind = which(!is.na(step))
  for(j in 1:2) allprobs[ind,j] = dgamma2(step[ind], mu[j], sigma[j])
  -forward(delta, Gamma, allprobs)
}

par = list(logitGamma = c(-2,-2), 
           logMu = log(c(0.3, 2.5)), 
           logSigma = log(c(0.3, 0.5)))
           
obj = MakeADFun(nll, par)
opt = nlminb(obj$par, obj$fn, obj$gr)

mod = obj$report()

pres = pseudo_res(step, "gamma2", list(mean = mod$mu, sd = mod$sigma),
                  mod = mod)
                  
plot(pres)
plot(pres, hist = TRUE)
}
