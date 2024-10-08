% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tpm_g.R
\name{tpm_g}
\alias{tpm_g}
\title{Build all transition probability matrices of an inhomogeneous HMM}
\usage{
tpm_g(Z, beta, byrow = FALSE, ad = FALSE, report = TRUE)
}
\arguments{
\item{Z}{Covariate design matrix with or without intercept column, i.e. of dimension c(n, p) or c(n, p+1).
If Z has only p columns, an intercept column of ones will be added automatically.}

\item{beta}{Matrix of coefficients for the off-diagonal elements of the transition probability matrix.
Needs to be of dimension c(N*(N-1), p+1), where the first column contains the intercepts.}

\item{byrow}{Logical that indicates if each transition probability matrix should be filled by row. 
Defaults to FALSE, but should be set to TRUE if one wants to work with a matrix of beta parameters returned by popular HMM packages like moveHMM, momentuHMM, or hmmTMB.}

\item{ad}{Logical, indicating whether automatic differentiation with RTMB should be used. Defaults to FALSE.}

\item{report}{Logical, indicating whether the coefficient matrix beta should be reported from the fitted model. Defaults to TRUE, but only works if ad = TRUE.}
}
\value{
Array of transition probability matrices of dimension c(N,N,n)
}
\description{
In an HMM, we model the influence of covariates on the state process, by linking them to the transition probabiltiy matrix. 
Most commonly, this is done by specifying a linear predictor
\deqn{ \eta_{ij}^{(t)} = \beta^{(ij)}_0 + \beta^{(ij)}_1 z_{t1} + \dots + \beta^{(ij)}_p z_{tp} }
for each off-diagonal element (\eqn{i \neq j}) and then applying the inverse multinomial logistic link to each row.
This function efficiently calculates all transition probabilty matrices for a given design matrix \eqn{Z} and parameter matrix beta.
}
\examples{
n = 1000
Z = matrix(runif(n*2), ncol = 2)
beta = matrix(c(-1, 1, 2, -2, 1, -2), nrow = 2, byrow = TRUE)
Gamma = tpm_g(Z, beta)
}
