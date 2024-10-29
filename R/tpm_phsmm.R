#' Build all transition probability matrices of an periodic-HSMM-approximating HMM
#'
#' Hidden semi-Markov models (HSMMs) are a flexible extension of HMMs. For direct numerical maximum likelhood estimation, HSMMs can be represented as HMMs on an enlarged state space (of size \eqn{M}) and with structured transition probabilities.
#' This function computes the transition matrices of a periodically inhomogeneos HSMMs.
#'
#' @param omega embedded transition probability matrix
#' 
#' Either a matrix of dimension c(N,N) for homogeneous conditional transition probabilities, or an array of dimension c(N,N,L) for inhomogeneous conditional transition probabilities.
#' @param dm state dwell-time distributions arranged in a list of length(N)
#' 
#' Each list element needs to be a matrix of dimension c(L, N_i), where each row t is the (approximate) probability mass function of state i at time t.
#' @param eps rounding value: If an entry of the transition probabily matrix is smaller, than it is rounded to zero.
#'
#' @return array of dimension c(N,N,L), containing the extended-state-space transition probability matrices of the approximating HMM for each time point of the cycle.
#' @export
#'
#' @examples
#' N = 3
#' L = 24
#' # time-varying mean dwell times
#' Lambda = exp(matrix(rnorm(L*N, 2, 0.5), nrow = L))
#' sizes = c(25, 25, 25) # approximating chain with 75 states
#' # state dwell-time distributions
#' dm = list()
#' for(i in 1:3){
#'   dmi = matrix(nrow = L, ncol = sizes[i])
#'   for(t in 1:L){
#'     dmi[t,] = dpois(1:sizes[i]-1, Lambda[t,i])
#'   }
#'   dm[[i]] = dmi
#' }
#' 
#' ## homogeneous conditional transition probabilites
#' # diagonal elements are zero, rowsums are one
#' omega = matrix(c(0,0.5,0.5,0.2,0,0.8,0.7,0.3,0), nrow = N, byrow = TRUE)
#' 
#' # calculating extended-state-space t.p.m.s
#' Gamma = tpm_phsmm(omega, dm)
#' 
#' ## inhomogeneous conditional transition probabilites
#' # omega can be an array
#' omega = array(rep(omega,L), dim = c(N,N,L))
#' omega[1,,4] = c(0, 0.2, 0.8) # small change for inhomogeneity
#' 
#' # calculating extended-state-space t.p.m.s
#' Gamma = tpm_phsmm(omega, dm)
tpm_phsmm2 = function(omega, dm, eps = 1e-10){
  # dm list over states: entries matrices of dim c(L, N_i)
  L = nrow(dm[[1]]) # length of one cycle
  N = length(dm) # number of states
  # allowing for inhomogeneous cond. transition probs
  if(is.matrix(omega)){
    omega = array(rep(omega,L), dim = c(N,N,L))
  }
  mv = sapply(dm, ncol) # lengths of the pmf vectors
  M = sum(mv)
  G_all = array(dim = c(M,M,L))
  Fm = vector("list") # computing all cdfs
  for(i in 1:N){
    Fm[[i]] = cbind(0,t(apply(dm[[i]][,-mv[i]], 1, cumsum)))
  }
  for (t in 1:L){ # loop over time
    G = matrix(0,0,sum(mv)) # constructing an empty matrix with number of columns equal to final dimension (sum_i m_i)
    for (i in 1:N){ # loop over states
      Ni = mv[i]
      ci = numeric(Ni) # length of pmf for state aggregate i
      for (k in 0:(Ni-1)){
        l = ifelse(t==k, L, (t-k)%%L)
        l = ifelse(l==0, L, l)
        ci[k+1] = ifelse(abs(1-Fm[[i]][l,k+1]) > eps,
                         dm[[i]][l,k+1]/(1-Fm[[i]][l,k+1]), 1)
      }
      cim = ifelse(1-ci > 0, 1-ci, 0)
      Gi = matrix(0, Ni, 0)
      for (j in 1:N){
        if(i==j) {
          if(Ni==1){ 
            Gi = cbind(Gi, c(rep(0,mv[[j]]-1), cim))
          } else{ 
            Gi = cbind(Gi, rbind(cbind(rep(0,Ni-1),diag(cim[-Ni],Ni-1,Ni-1)),
                                c(rep(0,Ni-1),cim[[Ni]])))}
        } else { if(Ni==1)
        { Gi = cbind(Gi, matrix(c(omega[[i,j,t]]*ci, rep(0,mv[[j]]-1)),1))} else
        { Gi = cbind(Gi, cbind(omega[[i,j,t]]*ci, matrix(0,mv[[i]],mv[[j]]-1)))}
        }
      }
      G = rbind(G,Gi)
    }
    G_all[,,t] = G
  }
  G_all
}