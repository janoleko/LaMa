% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generator.R
\name{generator}
\alias{generator}
\title{Build the generator matrix of a continuous-time Markov chain}
\usage{
generator(param, byrow = FALSE, report = TRUE)
}
\arguments{
\item{param}{Unconstraint parameter vector of length N*(N-1) where N is the number of states of the Markov chain}

\item{byrow}{Logical that indicates if the transition probability matrix should be filled by row. 
Defaults to FALSE, but should be set to TRUE if one wants to work with a matrix of beta parameters returned by popular HMM packages like moveHMM, momentuHMM, or hmmTMB.}

\item{report}{Logical, indicating whether the generator matrix Q should be reported from the fitted model. Defaults to TRUE, but only works if ad = TRUE.}
}
\value{
Infinitesimal generator matrix of dimension c(N,N)
}
\description{
This function builds the infinitesimal generator matrix for a continuous-time Markov chain from an unconstraint parameter vector.\cr

Compatible with automatic differentiation by RTMB
}
\examples{
# 2 states: 2 free off-diagonal elements
generator(rep(-1, 2))
# 3 states: 6 free off-diagonal elements
generator(rep(-2, 6))
}
