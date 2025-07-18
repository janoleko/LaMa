% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tpm_functions.R
\name{tpm_g}
\alias{tpm_g}
\title{Build all transition probability matrices of an inhomogeneous HMM}
\usage{
tpm_g(Z, beta, byrow = FALSE, ad = NULL, report = TRUE)
}
\arguments{
\item{Z}{covariate design matrix with or without intercept column, i.e. of dimension c(n, p) or c(n, p+1)

If \code{Z} has only p columns, an intercept column of ones will be added automatically.}

\item{beta}{matrix of coefficients for the off-diagonal elements of the transition probability matrix

Needs to be of dimension c(N*(N-1), p+1), where the first column contains the intercepts.}

\item{byrow}{logical indicating if each transition probability matrix should be filled by row
 
Defaults to \code{FALSE}, but should be set to \code{TRUE} if one wants to work with a matrix of beta parameters returned by popular HMM packages like \code{moveHMM}, \code{momentuHMM}, or \code{hmmTMB}.}

\item{ad}{optional logical, indicating whether automatic differentiation with \code{RTMB} should be used. By default, the function determines this itself.}

\item{report}{logical, indicating whether the coefficient matrix \code{beta} should be reported from the fitted model. Defaults to \code{TRUE}, but only works if \code{ad = TRUE}.}
}
\value{
array of transition probability matrices of dimension c(N,N,n)
}
\description{
In an HMM, we often model the influence of covariates on the state process by linking them to the transition probabiltiy matrix. 
Most commonly, this is done by specifying a linear predictor
\deqn{ \eta_{ij}^{(t)} = \beta^{(ij)}_0 + \beta^{(ij)}_1 z_{t1} + \dots + \beta^{(ij)}_p z_{tp} }
for each off-diagonal element (\eqn{i \neq j}) of the transition probability matrix and then applying the inverse multinomial logistic link (also known as softmax) to each row.
This function efficiently calculates all transition probabilty matrices for a given design matrix \code{Z} and parameter matrix \code{beta}.
}
\examples{
Z = matrix(runif(200), ncol = 2)
beta = matrix(c(-1, 1, 2, -2, 1, -2), nrow = 2, byrow = TRUE)
Gamma = tpm_g(Z, beta)
}
\seealso{
Other transition probability matrix functions: 
\code{\link{generator}()},
\code{\link{tpm}()},
\code{\link{tpm_cont}()},
\code{\link{tpm_emb}()},
\code{\link{tpm_emb_g}()},
\code{\link{tpm_g2}()},
\code{\link{tpm_p}()}
}
\concept{transition probability matrix functions}
