#' Compute the periodically stationary distribution of a periodically inhomogeneous Markov chain
#'
#' If the transition probability matrix of an inhomogeneous Markov chain varies only periodically (with period length \eqn{L}), it converges to a so-called periodically stationary distribution. 
#' This happens, because the thinned Markov chain, which has a full cycle as each time step, has homogeneous transition probability matrix \cr\cr
#' \eqn{\Gamma_t = \Gamma^{(t)} \Gamma^{(t+1)} \dots \Gamma^{(t+L-1)}} for all \eqn{t = 1, \dots, L}. \cr \cr
#' The stationary distribution for time \eqn{t} satifies \eqn{\delta^{(t)} \Gamma_t = \delta^{(t)}}. \cr
#' This function calculates the periodically stationary distribution.
#' 
#' @param Gamma Array of transition probability matrices of dimension c(N,N,L). 
#' @param t Integer index of the time point in the cycle, for which to calculate the stationary distribution
#' If t is not provided, the function calculates all stationary distributions for each time point in the cycle.
#' @param tol The tolerance for detecting linear dependencies in the columns of the thinned transition matrix. The default is .Machine$double.eps.
#'
#' @return Either the periodically stationary distribution at time t or all periodically stationary distributions.
#' @export
#'
#' @examples
# setting parameters for trigonometric link
#' beta = matrix(c(-1, -2, 2, -1, 2, -4), nrow = 2, byrow = TRUE)
#' # building trigonometric design matrix
#' Z = cbind(1,trigBasisExp(1:24, 24, 1))
#' # calculating all 24 linear predictor vectors
#' Eta = Z%*%t(beta)
#' # building all 24 t.p.m.s
#' Gamma = array(dim = c(2,2,24))
#' for(t in 1:24){
#'   Gamma[,,t] = tpm(Eta[t,])
#' }
#' # Periodically stationary distribution for specific time point
#' delta4 = stationary_p(Gamma, 4)
#'
#' # All periodically stationary distributions
#' Delta = stationary_p(Gamma)
stationary_p = function(Gamma, t = NULL, tol = .Machine$double.eps){
  if(is.null(t)){
    N = dim(Gamma)[2]
    L = dim(Gamma)[3]
    Delta = matrix(nrow = L, ncol = N)
    GammaT = tpm_thinned(Gamma, 1)
    Delta[1,] = stationary(GammaT, tol)
    for(t in 2:L){
      Delta[t,] = Delta[t-1,]%*%Gamma[,,t-1]
    }
    return(Delta)
  } else{
    GammaT = tpm_thinned(Gamma, t)
    delta = stationary(GammaT, tol)
    return(delta)
  }
}
