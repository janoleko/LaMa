% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/decoding_functions.R
\name{stateprobs}
\alias{stateprobs}
\title{Calculate conditional local state probabilities for homogeneous HMMs}
\usage{
stateprobs(delta, Gamma, allprobs, trackID = NULL, mod = NULL)
}
\arguments{
\item{delta}{initial or stationary distribution of length N, or matrix of dimension c(k,N) for k independent tracks, if \code{trackID} is provided}

\item{Gamma}{transition probability matrix of dimension c(N,N), or array of k transition probability matrices of dimension c(N,N,k), if \code{trackID} is provided}

\item{allprobs}{matrix of state-dependent probabilities/ density values of dimension c(n, N)}

\item{trackID}{optional vector of length n containing IDs

If provided, the total log-likelihood will be the sum of each track's likelihood contribution.
In this case, \code{Gamma} can be a matrix, leading to the same transition probabilities for each track, or an array of dimension c(N,N,k), with one (homogeneous) transition probability matrix for each track.
Furthermore, instead of a single vector \code{delta} corresponding to the initial distribution, a \code{delta} matrix of initial distributions, of dimension c(k,N), can be provided, such that each track starts with it's own initial distribution.}

\item{mod}{optional model object containing initial distribution \code{delta}, transition probability matrix \code{Gamma}, matrix of state-dependent probabilities \code{allprobs}, and potentially a \code{trackID} variable

If you are using automatic differentiation either with \code{RTMB::MakeADFun} or \code{\link{qreml}} and include \code{\link{forward}} in your likelihood function, the objects needed for state decoding are automatically reported after model fitting.
Hence, you can pass the model object obtained from running \code{RTMB::report()} or from \code{\link{qreml}} directly to this function.}
}
\value{
matrix of conditional state probabilities of dimension c(n,N)
}
\description{
Computes
\deqn{\Pr(S_t = j \mid X_1, ..., X_T)}
for homogeneous HMMs
}
\examples{
Gamma = tpm(c(-1,-2))
delta = stationary(Gamma)
allprobs = matrix(runif(10), nrow = 10, ncol = 2)

probs = stateprobs(delta, Gamma, allprobs)
}
\seealso{
Other decoding functions: 
\code{\link{stateprobs_g}()},
\code{\link{stateprobs_p}()},
\code{\link{viterbi}()},
\code{\link{viterbi_g}()},
\code{\link{viterbi_p}()}
}
\concept{decoding functions}
