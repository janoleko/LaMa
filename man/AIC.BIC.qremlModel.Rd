% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/qreml_functions.R
\name{AIC.BIC.qremlModel}
\alias{AIC.BIC.qremlModel}
\alias{AIC.qremlModel}
\alias{BIC.qremlModel}
\title{Akaike's Information Criterion for models returned by \code{\link{qreml}}}
\usage{
\method{AIC}{qremlModel}(object, k = 2)

\method{BIC}{qremlModel}(object, ...)
}
\arguments{
\item{object}{a fitted model object returned by \code{\link{qreml}}}

\item{k}{numeric, the penalty per parameter to be used; the default \code{k = 2} is the classical AIC.}

\item{...}{can include \code{nObs} number of observations, only used for BIC. If not provided, the function tries to extract it from the model object.}
}
\value{
numeric value with the corresponding AIC (or BIC, or ..., depending on k)
}
\description{
Akaike's Information Criterion for models returned by \code{\link{qreml}}
}
\examples{
## no examples
}
