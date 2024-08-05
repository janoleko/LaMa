#' Build the transition probability matrix from unconstraint parameter vector
#'
#' @description
#' This function builds the transition probability matrix from an unconstraint parameter vector. 
#' For each row of the matrix, the inverse multinomial logistic link is applied.
#'
#' @param param Unconstraint parameter vector of length N*(N-1) where N is the number of states of the Markov chain
#' @param byrow Logical that indicates if the transition probability matrix should be filled by row. 
#' Defaults to FALSE, but should be set to TRUE if one wants to work with a matrix of beta parameters returned by popular HMM packages like moveHMM, momentuHMM, or hmmTMB.
#'
#' @return Transition probability matrix of dimension c(N,N)
#' @export
#'
#' @examples
#' # 2 states: 2 free off-diagonal elements
#' param1 = rep(-1, 2)
#' Gamma1 = tpm(param1)
#' 
#' # 3 states: 6 free off-diagonal elements
#' param2 = rep(-2, 6)
#' Gamma2 = tpm(param2)
tpm = function(param, byrow = FALSE) {
  K = length(param)
  # for N > 1: N*(N-1) is bijective with solution
  N = as.integer(0.5 + sqrt(0.25 + K), 0)
  
  Gamma = diag(N)
  Gamma[!Gamma] = exp(param[1:(N*(N-1))])
  if(byrow == TRUE){
    Gamma = t(Gamma)
  }
  Gamma = Gamma / rowSums(Gamma)
  Gamma
}
