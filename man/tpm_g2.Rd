% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tpm_functions.R
\name{tpm_g2}
\alias{tpm_g2}
\title{Build all transition probability matrices of an inhomogeneous HMM}
\usage{
tpm_g2(Z, beta, byrow = FALSE, ad = NULL, report = TRUE, ref = NULL)
}
\arguments{
\item{Z}{covariate design matrix with or without intercept column, i.e. of dimension c(n, p) or c(n, p+1)

If \code{Z} has only p columns, an intercept column of ones will be added automatically.

Can also be a list of N*(N-1) design matrices with different number of columns but the same number of rows. In that case, no intercept column will be added.}

\item{beta}{matrix of coefficients for the off-diagonal elements of the transition probability matrix

Needs to be of dimension c(N*(N-1), p+1), where the first column contains the intercepts.

If \code{Z} is a list, \code{beta} can also be a list of length N*(N-1) with each entry being a vector or a (long) matrix of coefficients, each matching the dimension of the corresponding entry in \code{Z}.}

\item{byrow}{logical indicating if each transition probability matrix should be filled by row
 
Defaults to \code{FALSE}, but should be set to \code{TRUE} if one wants to work with a matrix of beta parameters returned by popular HMM packages like \code{moveHMM}, \code{momentuHMM}, or \code{hmmTMB}.}

\item{ad}{optional logical, indicating whether automatic differentiation with \code{RTMB} should be used. By default, the function determines this itself.}

\item{report}{logical, indicating whether the coefficient matrix \code{beta} should be reported from the fitted model. Defaults to \code{TRUE}, but only works if \code{ad = TRUE}.}

\item{ref}{optional vector of length N with the reference state indices for each column of the transition probability matrix. Each row in the transition matrix corresponds to a multinomial regression, hence one state needs to be the reference category. Defaults to off-diagonal elements (\code{ref = 1:N}).}
}
\value{
array of transition probability matrices of dimension c(N,N,n)
}
\description{
In an HMM, we often model the influence of covariates on the state process by linking them to the transition probabiltiy matrix. 
Most commonly, this is done by specifying linear predictors
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
\code{\link{tpm_g}()},
\code{\link{tpm_p}()}
}
\concept{transition probability matrix functions}
