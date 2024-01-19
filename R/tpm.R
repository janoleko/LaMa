#' Build the transition probability matrix from unconstraint parameter vector
#'
#' @description
#' This function builds the transition probability matrix from an unconstraint parameter vector. 
#' For each row of the matrix, the inverse multinomial logistic link is applied.
#'
#' @param param Unconstraint parameter vector of length N*(N-1) where N is the number of states of the Markov chain
#'
#' @return Transition probability matrix of dimension c(N,N)
#' @export
#'
#' @examples
#' param = rep(-2,6)
#' Gamma = tpm(param)
tpm = function(param){
  K = length(param)
  # for N > 1: N*(N-1) is bijective with solution
  N = as.integer(0.5 + sqrt(0.25+K), 0)
  
  Gamma = diag(N)
  Gamma[!Gamma] = exp(param[1:(N*(N-1))])
  Gamma = Gamma / rowSums(Gamma)
  Gamma
}
