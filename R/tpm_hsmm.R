#' Build the transition probability matrix of an HSMM-approximating HMM
#'
#' @description
#' Hidden semi-Markov models (HSMMs) are a flexible extension of HMMs. 
#' For direct numerical maximum likelhood estimation, HSMMs can be represented as HMMs on an enlarged state space (of size \eqn{M}) and with structured transition probabilities.
#' This function computes the transition matrix of an HSMM.
#'
#' @param omega Embedded transition probability matrix of dimension c(N,N)
#' @param dm State dwell-time distributions arranged in a list of length(N). Each list element needs to be a vector of length N_i, where N_i is the state aggregate size.
#' @param eps Rounding value: If an entry of the transition probabily matrix is smaller, than it is rounded to zero.
#'
#' @return The extended-state-space transition probability matrix of the approximating HMM
#' @export
#'
#' @examples
#' # building the t.p.m. of the embedded Markov chain
#' omega = matrix(c(0,1,1,0), nrow = 2, byrow = TRUE)
#' # defining state aggregate sizes
#' sizes = c(20, 30)
#' # defining state dwell-time distributions
#' lambda = c(5, 11)
#' dm = list(dpois(1:sizes[1]-1, lambda[1]), dpois(1:sizes[2]-1, lambda[2]))
#' # calculating extended-state-space t.p.m.
#' Gamma = tpm_hsmm(omega, dm)
tpm_hsmm = function(omega,dm,eps=1e-10){
  mv = sapply(dm,length)
  m = length(mv)
  G = matrix(0,0,sum(mv))
  for (i in 1:m){
    mi = mv[[i]]
    F = cumsum(c(0,dm[[i]][-mi]))
    ci = ifelse(abs(1-F)>eps,dm[[i]]/(1-F),1)
    cim = ifelse(1-ci>0,1-ci,0)
    Gi = matrix(0,mi,0)
    for (j in 1:m){
      if(i==j) {
        if(mi==1){ 
          Gi = cbind(Gi,c(rep(0,mv[[j]]-1),cim))
        } else{ 
          Gi = cbind(Gi,rbind(cbind(rep(0,mi-1),diag(cim[-mi],mi-1,mi-1)),
                              c(rep(0,mi-1),cim[[mi]])))}
      } else   { if(mi==1)
      { Gi = cbind(Gi,matrix(c(omega[[i,j]]*ci,rep(0,mv[[j]]-1)),1))} else
      { Gi = cbind(Gi,cbind(omega[[i,j]]*ci,matrix(0,mv[[i]],mv[[j]]-1)))}
      }
    }
    G = rbind(G,Gi)
  }
  G 
}