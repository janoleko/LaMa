% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stationary_functions.R
\name{stationary_cont}
\alias{stationary_cont}
\title{Compute the stationary distribution of a continuous-time Markov chain}
\usage{
stationary_cont(Q)
}
\arguments{
\item{Q}{infinitesimal generator matrix of dimension \code{c(N,N)} or array of such matrices of dimension \code{c(N,N,nTracks)} if the stationary distribution should be computed for several matrices at once}
}
\value{
either a single stationary distribution of the continuous-time Markov chain (vector of length \code{N}) or a matrix of stationary distributions of dimension \code{c(nTracks,N)} with one stationary distribution in each row
}
\description{
A well-behaved continuous-time Markov chain converges to a unique stationary distribution, here called \eqn{\pi}.
This distribution satisfies
\deqn{\pi Q = 0,} subject to \eqn{\sum_{j=1}^N \pi_j = 1},
where \eqn{Q} is the infinitesimal generator of the Markov chain.
This function solves the linear system of equations above for a given generator matrix.
}
\examples{
# single matrix
Q = generator(c(-2,-2))
Pi = stationary_cont(Q)
# multiple matrices
Q = array(Q, dim = c(2,2,10))
Pi = stationary_cont(Q)
}
\seealso{
\code{\link{generator}} to create a generator matrix

Other stationary distribution functions: 
\code{\link{stationary}()},
\code{\link{stationary_p}()}
}
\concept{stationary distribution functions}
