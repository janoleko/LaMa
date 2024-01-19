#' Compute the transition probability matrix of a thinned periodically inhomogeneous Markov chain.
#'
#' If the transition probability matrix of an inhomogeneous Markov chain varies only periodically (with period length \eqn{L}), it converges to a so-called periodically stationary distribution. 
#' This happens, because the thinned Markov chain, which has a full cycle as each time step, has homogeneous transition probability matrix \cr\cr
#' \eqn{\Gamma_t = \Gamma^{(t)} \Gamma^{(t+1)} \dots \Gamma^{(t+L-1)}} for all \eqn{t = 1, \dots, L}. \cr \cr
#' This function calculates the matrix above efficiently as a preliminery step to calculating the periodically stationary distribution.
#'
#' @param Gamma Array of transition probability matrices of dimension c(N,N,L).
#' @param t Integer index of the time point in the cycle, for which to calculate the thinned transition probility matrix
#'
#' @return Thinned transition probabilty matrix of dimension c(N,N)
#' @export
#'
#' @examples
#' # setting parameters for trigonometric link
#' beta = matrix(c(-1, -2, 2, -1, 2, -4), nrow = 2, byrow = TRUE)
#' # building trigonometric design matrix
#' Z = cbind(1,trigBasisExp(1:24, 24, 1))
#' # calculating all 24 linear predictor vectors
#' Eta = Z%*%t(beta)
#' # building all 24 t.p.m.s
#' Gamma = array(dim = c(2,2,24))
#' for(t in 1:24){
#'   Gamma[,,t] = tpm(Eta[t,], 2)
#' }
#' # calculating 
#' tpm_thinned(Gamma, 4)
tpm_thinned = function(Gamma, t){
  tpm_thinned_t_cpp(Gamma, t)
}



