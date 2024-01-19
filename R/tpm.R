#' Build the transition probability matrix from unconstraint parameter vector
#'
#' @description
#' This function builds the transition probability matrix from an unconstraint parameter vector. 
#' For each row of the matrix, the inverse multinomial logistic link is applied.
#'
#' @param param Unconstraint parameter vector of length N*(N-1)
#' @param N Number of hidden states
#'
#' @return Transition probability matrix of dimension c(N,N)
#' @export
#'
#' @examples
#' param = rep(-2,6)
#' Gamma = tpm(param, 3)
tpm = function(param, N){
  if(length(param)!= N*(N-1)){
    return(NA)
  }
  Gamma = diag(N)
  Gamma[!Gamma] = exp(param[1:(N*(N-1))])
  Gamma = Gamma / rowSums(Gamma)
  Gamma
}