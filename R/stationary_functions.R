#' Compute the stationary distribution of a homogeneous Markov chain
#'
#' @description
#' A homogeneous, finite state Markov chain that is irreducible and aperiodic converges to a unique stationary distribution, here called \eqn{\delta}.
#' As it is stationary, this distribution satisfies
#' \deqn{\delta \Gamma = \delta,} subject to \eqn{\sum_{j=1}^N \delta_j = 1},
#' where \eqn{\Gamma} is the transition probability matrix. 
#' This function solves the linear system of equations above.
#
#' @param Gamma transition probability matrix of dimension c(N,N)
#'
#' @return stationary distribution of the Markov chain with the given transition probability matrix
#' @export
#' @import RTMB
#'
#' @examples
#' Gamma = tpm(c(rep(-2,3), rep(-3,3)))
#' delta = stationary(Gamma)
stationary = function(Gamma){
  "[<-" <- ADoverload("[<-") # overloading assignment operators, currently necessary
  "c" <- ADoverload("c")
  "diag<-" <- ADoverload("diag<-")
  
  N = dim(Gamma)[1]
  delta = RTMB::solve(t(diag(N)-Gamma+1), rep(1,N))
  names(delta) = paste("state", 1:N)
  delta
}


#' Compute the periodically stationary distribution of a periodically inhomogeneous Markov chain
#'
#' @description
#' If the transition probability matrix of an inhomogeneous Markov chain varies only periodically (with period length \eqn{L}), it converges to a so-called periodically stationary distribution. 
#' This happens, because the thinned Markov chain, which has a full cycle as each time step, has homogeneous transition probability matrix
#' \deqn{\Gamma_t = \Gamma^{(t)} \Gamma^{(t+1)} \dots \Gamma^{(t+L-1)}} for all \eqn{t = 1, \dots, L.}
#' The stationary distribution for time \eqn{t} satifies \eqn{\delta^{(t)} \Gamma_t = \delta^{(t)}}.
#' 
#' This function calculates said periodically stationary distribution.
#' 
#' @param Gamma array of transition probability matrices of dimension c(N,N,L)
#' @param t integer index of the time point in the cycle, for which to calculate the stationary distribution
#' 
#' If t is not provided, the function calculates all stationary distributions for each time point in the cycle.
#' @param ad optional logical, indicating whether automatic differentiation with \code{RTMB} should be used. By default, the function determines this itself.
#'
#' @return either the periodically stationary distribution at time t or all periodically stationary distributions.
#' @export
#' @import RTMB
#'
#' @examples
# setting parameters for trigonometric link
#' L = 24
#' beta = matrix(c(-1, 2, -1, -2, 1, -1), nrow = 2, byrow = TRUE)
#' Gamma = tpm_p(1:L, L, beta, degree = 1)
#' # Periodically stationary distribution for specific time point
#' delta = stationary_p(Gamma, 4)
#'
#' # All periodically stationary distributions
#' Delta = stationary_p(Gamma)
stationary_p = function(Gamma, t = NULL, ad = NULL){
  
  N = dim(Gamma)[2]
  L = dim(Gamma)[3]
  
  # if ad is not explicitly provided, check if delta is an advector
  if(is.null(ad)){
    # check if delta has any of the allowed classes
    if(!any(class(Gamma) %in% c("advector", "array"))){
      stop("Gamma needs to be either an array or advector.")
    }
    
    # if delta is advector, run ad version of the function
    ad = inherits(Gamma, "advector")
  }
  
  if(!ad) {
    if(is.null(t)){
      Delta = matrix(nrow = L, ncol = N)
      GammaT = tpm_thinned(Gamma, 1)
      Delta[1,] = stationary(GammaT)
      for(t in 2:L){
        Delta[t,] = Delta[t-1,]%*%Gamma[,,t-1]
      }
      colnames(Delta) = paste("state", 1:N)
      return(Delta)
    } else{
      GammaT = tpm_thinned(Gamma, t)
      delta = stationary(GammaT)
      names(delta) = paste("state", 1:length(delta))
      return(delta)
    }
  } else if(ad) {
    
    "[<-" <- ADoverload("[<-") # overloading assignment operators, currently necessary
    "c" <- ADoverload("c")
    "diag<-" <- ADoverload("diag<-")
    
    if(is.null(t)) {
      Delta = matrix(NaN, nrow = L, ncol = N)
      
      GammaT = Gamma[,,1]
      for(k in 2:L) GammaT = GammaT %*% Gamma[,,k]
      
      Delta[1,] = stationary(GammaT)
      
      for(t in 2:L){
        Delta[t,] = c(t(Delta[t-1,]) %*% Gamma[,,t-1])
      }
      colnames(Delta) = paste("state", 1:N)
      return(Delta)
    } else{
      GammaT = Gamma[,,t]
      if(t < L){
        for(k in (t+1):L) GammaT = GammaT %*% Gamma[,,k]
        if(t > 1){
          for(k in 1:(t-1)) GammaT = GammaT %*% Gamma[,,k]
        }
      } else if(t == L){
        for(k in 1:(L-1)) GammaT = GammaT %*% Gamma[,,k]
      }
      delta = stationary(GammaT)
      return(delta)
    }
  }
}


#' Compute the stationary distribution of a continuous-time Markov chain
#'
#' @description
#' A well-behaved continuous-time Markov chain converges to a unique stationary distribution, here called \eqn{\pi}.
#' This distribution satisfies
#' \deqn{\pi Q = 0,} subject to \eqn{\sum_{j=1}^N \pi_j = 1},
#' where \eqn{Q} is the infinitesimal generator of the Markov chain.
#' This function solves the linear system of equations above.
#
#' @param Q infinitesimal generator matrix of dimension c(N,N)
#'
#' @return stationary distribution of the continuous-time Markov chain with given generator matrix
#' @export
#' @import RTMB
#'
#' @examples
#' Q = generator(c(-2,-2))
#' Pi = stationary_cont(Q)
stationary_cont = function(Q){
  "[<-" <- ADoverload("[<-") # overloading assignment operators, currently necessary
  "c" <- ADoverload("c")
  "diag<-" <- ADoverload("diag<-")
  
  N = dim(Q)[1]
  Pi = RTMB::solve(t(Q + 1), rep(1,N))
  names(Pi) = paste("state", 1:N)
  Pi
}


#' Sparse version of \code{\link{stationary}}
#'
#' @description
#' This is function computes the stationary distribution of a Markov chain with a given \strong{sparse} transition probability matrix.
#' Compatible with automatic differentiation by \code{RTMB}
#
#' @param Gamma sparse transition probability matrix of dimension c(N,N)
#'
#' @return stationary distribution of the Markov chain with the given transition probability matrix
#' @export
#' @import RTMB
#'
#' @examples
#' ## HSMM example (here the approximating tpm is sparse)
#' # building the t.p.m. of the embedded Markov chain
#' omega = matrix(c(0,1,1,0), nrow = 2, byrow = TRUE)
#' # defining state aggregate sizes
#' sizes = c(20, 30)
#' # defining state dwell-time distributions
#' lambda = c(5, 11)
#' dm = list(dpois(1:sizes[1]-1, lambda[1]), dpois(1:sizes[2]-1, lambda[2]))
#' # calculating extended-state-space t.p.m.
#' Gamma = tpm_hsmm(omega, dm)
#' delta = stationary_sparse(Gamma)
stationary_sparse = function(Gamma) 
{
  "[<-" <- ADoverload("[<-")
  "c" <- ADoverload("c")
  "diag<-" <- ADoverload("diag<-")
  N = dim(Gamma)[1]
  
  if(methods::is(Gamma, "sparseMatrix")) {
    Gamma = as.matrix(Gamma)
  }
  delta = Matrix::solve(t(diag(N) - Gamma + 1), rep(1, N))
  
  delta
}

#' Sparse version of \code{\link{stationary_p}}
#'
#' @description
#' This is function computes the periodically stationary distribution of a Markov chain given a list of L \strong{sparse} transition probability matrices.
#' Compatible with automatic differentiation by \code{RTMB}
#' 
#' @param Gamma sist of length L containing sparse transition probability matrices for one cycle.
#' @param t integer index of the time point in the cycle, for which to calculate the stationary distribution
#' If t is not provided, the function calculates all stationary distributions for each time point in the cycle.
#'
#' @return either the periodically stationary distribution at time t or all periodically stationary distributions.
#' @export
#' @import RTMB
#'
#' @examples
#' ## periodic HSMM example (here the approximating tpm is sparse)
#' N = 2 # number of states
#' L = 24 # cycle length
#' # time-varying mean dwell times
#' Z = trigBasisExp(1:L) # trigonometric basis functions design matrix
#' beta = matrix(c(2, 2, 0.1, -0.1, -0.2, 0.2), nrow = 2)
#' Lambda = exp(cbind(1, Z) %*% t(beta))
#' sizes = c(20, 20) # approximating chain with 40 states
#' # state dwell-time distributions
#' dm = lapply(1:N, function(i) sapply(1:sizes[i]-1, dpois, lambda = Lambda[,i]))
#' omega = matrix(c(0,1,1,0), nrow = N, byrow = TRUE) # embedded t.p.m.
#' 
#' # calculating extended-state-space t.p.m.s
#' Gamma = tpm_phsmm(omega, dm)
#' # Periodically stationary distribution for specific time point
#' delta = stationary_p_sparse(Gamma, 4)
#'
#' # All periodically stationary distributions
#' Delta = stationary_p_sparse(Gamma)
stationary_p_sparse = function(Gamma, t = NULL){
  # overloading assignment operators, currently necessary
  "[<-" <- ADoverload("[<-")
  "c" <- ADoverload("c")
  "diag<-" <- ADoverload("diag<-")
  
  L = length(Gamma) # cycle length
  N = dim(Gamma[[1]])[1] # tpm dimension
  
  if(is.null(t)) {
    Delta = matrix(NaN, nrow = L, ncol = N)
    
    GammaT = Gamma[[1]]
    for(k in 2:L) GammaT = GammaT %*% Gamma[[k]]
    
    Delta[1,] = stationary_sparse(GammaT)
    
    for(t in 2:L){
      Delta[t,] = as.numeric(Delta[t-1,, drop = FALSE] %*% Gamma[[t-1]])
    }
    return(Delta)
  } else{
    GammaT = Gamma[[t]]
    if(t < L){
      for(k in (t+1):L) GammaT = GammaT %*% Gamma[[k]]
      if(t > 1){
        for(k in 1:(t-1)) GammaT = GammaT %*% Gamma[[k]]
      }
    } else if(t == L){
      for(k in 1:(L-1)) GammaT = GammaT %*% Gamma[[k]]
    }
    delta = stationary_sparse(GammaT)
    return(delta)
  }
}