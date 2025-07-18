
# Functions for HMMs ------------------------------------------------------


#' Build the transition probability matrix from unconstrained parameter vector
#'
#' @description
#' Markov chains are parametrised in terms of a transition probability matrix \eqn{\Gamma}, for which each row contains a conditional probability distribution of the next state given the current state.
#' Hence, each row has entries between 0 and 1 that need to sum to one. 
#' 
#' For numerical optimisation, we parametrise in terms of unconstrained parameters, thus this function computes said matrix from an unconstrained parameter vector via the inverse multinomial logistic link (also known as softmax) applied to each row.
#'
#' @family transition probability matrix functions
#'
#' @param param unconstrained parameter vector of length N*(N-1) where N is the number of states of the Markov chain
#' @param byrow logical indicating if the transition probability matrix should be filled by row
#' 
#' Defaults to \code{FALSE}, but should be set to \code{TRUE} if one wants to work with a matrix of beta parameters returned by popular HMM packages like \code{moveHMM}, \code{momentuHMM}, or \code{hmmTMB}.
#'
#' @return Transition probability matrix of dimension c(N,N)
#' @export
#' @import RTMB
#'
#' @examples
#' # 2 states: 2 free off-diagonal elements
#' par1 = rep(-1, 2)
#' Gamma1 = tpm(par1)
#' 
#' # 3 states: 6 free off-diagonal elements
#' par2 = rep(-2, 6)
#' Gamma2 = tpm(par2)
tpm = function(param, byrow = FALSE) {
  
  "[<-" <- RTMB::ADoverload("[<-") # overloading assignment operators, currently necessary
  "c" <- ADoverload("c")
  "diag<-" <- ADoverload("diag<-")
  
  K = length(param)
  # for N > 1: N*(N-1) is bijective with solution
  N = as.integer(0.5 + sqrt(0.25 + K), 0)
  
  Gamma = diag(N)
  Gamma[!Gamma] = exp(param[1:(N*(N-1))])
  
  if(byrow) Gamma = t(Gamma)
  
  Gamma = Gamma / rowSums(Gamma)
  
  # naming
  statenames <- paste0("S", 1:N)
  rownames(Gamma) <- statenames
  colnames(Gamma) <- statenames
  
  Gamma
}


#' Build all transition probability matrices of an inhomogeneous HMM
#' 
#' In an HMM, we often model the influence of covariates on the state process by linking them to the transition probabiltiy matrix. 
#' Most commonly, this is done by specifying a linear predictor
#' \deqn{ \eta_{ij}^{(t)} = \beta^{(ij)}_0 + \beta^{(ij)}_1 z_{t1} + \dots + \beta^{(ij)}_p z_{tp} }
#' for each off-diagonal element (\eqn{i \neq j}) of the transition probability matrix and then applying the inverse multinomial logistic link (also known as softmax) to each row.
#' This function efficiently calculates all transition probabilty matrices for a given design matrix \code{Z} and parameter matrix \code{beta}.
#' 
#' @family transition probability matrix functions
#'
#' @param Z covariate design matrix with or without intercept column, i.e. of dimension c(n, p) or c(n, p+1)
#' 
#' If \code{Z} has only p columns, an intercept column of ones will be added automatically.
#' @param beta matrix of coefficients for the off-diagonal elements of the transition probability matrix
#' 
#' Needs to be of dimension c(N*(N-1), p+1), where the first column contains the intercepts.
#' @param byrow logical indicating if each transition probability matrix should be filled by row
#'  
#' Defaults to \code{FALSE}, but should be set to \code{TRUE} if one wants to work with a matrix of beta parameters returned by popular HMM packages like \code{moveHMM}, \code{momentuHMM}, or \code{hmmTMB}.
#' @param ad optional logical, indicating whether automatic differentiation with \code{RTMB} should be used. By default, the function determines this itself.
#' @param report logical, indicating whether the coefficient matrix \code{beta} should be reported from the fitted model. Defaults to \code{TRUE}, but only works if \code{ad = TRUE}.
#'
#' @return array of transition probability matrices of dimension c(N,N,n)
#' @export
#' @import RTMB
#' @importFrom stats na.omit
#'
#' @examples
#' Z = matrix(runif(200), ncol = 2)
#' beta = matrix(c(-1, 1, 2, -2, 1, -2), nrow = 2, byrow = TRUE)
#' Gamma = tpm_g(Z, beta)
tpm_g = function(Z, beta, byrow = FALSE, ad = NULL, report = TRUE){
  
  K = nrow(beta)
  p = ncol(beta) - 1
  N = as.integer(0.5 + sqrt(0.25 + K), 0)
  
  Z = as.matrix(Z)
  
  if(ncol(Z) == p){
    Z = cbind(1, Z) # adding intercept column
  } else if(ncol(Z) != p + 1){
    stop("The dimensions of Z and beta do not match.")
  }
  
  # report quantities for easy use later
  if(report) {
    #if(is.null(colnames(beta))){
      # Setting colnames for beta: Inherit colnames from Z
    colnames(beta) <- colnames(Z)
    #}
    if(is.null(rownames(beta))){
      # Setting rownames: depends on byrow
      names <- outer(paste0("S", 1:N, ">"), paste0("S", 1:N), FUN = paste0) # matrix
      diag(names) <- NA
      rownames(beta) <- na.omit(if (byrow) c(t(names)) else c(names))
    }
    RTMB::REPORT(beta)
  }
  
  # if ad is not explicitly provided, check if delta is an advector
  if(is.null(ad)){
    # check if delta has any of the allowed classes
    if(!any(class(beta) %in% c("advector", "numeric", "matrix", "array"))){
      stop("beta needs to be either a matrix or advector.")
    }
    
    # if delta is advector, run ad version of the function
    ad = inherits(beta, "advector")
  }
  
  if(!ad) {
    
    Gamma = tpm_g_cpp(Z, beta, N, byrow) # C++ version
    
  } else if(ad) {
    "[<-" <- ADoverload("[<-") # overloading assignment operators, currently necessary
    "c" <- ADoverload("c")
    "diag<-" <- ADoverload("diag<-")
    
    # if(report) {
    #   RTMB::REPORT(beta) # reporting coefficient matrix
    # }
    
    expEta = exp(Z %*% t(beta))
    Gamma = array(1, dim = c(N, N, nrow(expEta)))

    ## Loop over entries (stuff over time happens vectorised which speeds up the tape)
    col_ind <- 1
    for(i in seq_len(N)){ # loop over rows
      for(j in seq_len(N)){ # loop over columns
        if(j != i){ # only if non-diagonal
          if(byrow){
            Gamma[i, j, ] <- expEta[, col_ind]
          } else{
            Gamma[j, i, ] <- expEta[, col_ind]
          }
          # increase col_ind by one
          col_ind = col_ind + 1
        }
      }
    }
    
    # Normalise rows to sum to 1
    for(i in seq_len(N)){
      # transposing is necessary because Gamma[i,,] / colSums(Gamma[i,,]) does not work as expected
      Gamma[i, , ] <- t(t(Gamma[i, , ]) / rowSums(t(Gamma[i, , ])))
    }
  }
  
  # naming
  statenames <- paste0("S", 1:N)
  rownames(Gamma) <- statenames
  colnames(Gamma) <- statenames
  
  Gamma
}

#' Build all transition probability matrices of an inhomogeneous HMM
#' 
#' In an HMM, we often model the influence of covariates on the state process by linking them to the transition probabiltiy matrix. 
#' Most commonly, this is done by specifying linear predictors
#' \deqn{ \eta_{ij}^{(t)} = \beta^{(ij)}_0 + \beta^{(ij)}_1 z_{t1} + \dots + \beta^{(ij)}_p z_{tp} }
#' for each off-diagonal element (\eqn{i \neq j}) of the transition probability matrix and then applying the inverse multinomial logistic link (also known as softmax) to each row.
#' This function efficiently calculates all transition probabilty matrices for a given design matrix \code{Z} and parameter matrix \code{beta}.
#' 
#' @family transition probability matrix functions
#'
#' @param Z covariate design matrix with or without intercept column, i.e. of dimension c(n, p) or c(n, p+1)
#' 
#' If \code{Z} has only p columns, an intercept column of ones will be added automatically.
#' 
#' Can also be a list of N*(N-1) design matrices with different number of columns but the same number of rows. In that case, no intercept column will be added.
#' @param beta matrix of coefficients for the off-diagonal elements of the transition probability matrix
#' 
#' Needs to be of dimension c(N*(N-1), p+1), where the first column contains the intercepts.
#' 
#' If \code{Z} is a list, \code{beta} can also be a list of length N*(N-1) with each entry being a vector or a (long) matrix of coefficients, each matching the dimension of the corresponding entry in \code{Z}.
#' @param byrow logical indicating if each transition probability matrix should be filled by row
#'  
#' Defaults to \code{FALSE}, but should be set to \code{TRUE} if one wants to work with a matrix of beta parameters returned by popular HMM packages like \code{moveHMM}, \code{momentuHMM}, or \code{hmmTMB}.
#' @param ad optional logical, indicating whether automatic differentiation with \code{RTMB} should be used. By default, the function determines this itself.
#' @param report logical, indicating whether the coefficient matrix \code{beta} should be reported from the fitted model. Defaults to \code{TRUE}, but only works if \code{ad = TRUE}.
#' @param ref optional vector of length N with the reference state indices for each column of the transition probability matrix. Each row in the transition matrix corresponds to a multinomial regression, hence one state needs to be the reference category. Defaults to off-diagonal elements (\code{ref = 1:N}).
#'
#' @return array of transition probability matrices of dimension c(N,N,n)
#' @export
#' @import RTMB
#' @importFrom stats na.omit
#'
#' @examples
#' Z = matrix(runif(200), ncol = 2)
#' beta = matrix(c(-1, 1, 2, -2, 1, -2), nrow = 2, byrow = TRUE)
#' Gamma = tpm_g(Z, beta)
tpm_g2 <- function(Z, 
                   beta, 
                   byrow = FALSE, 
                   ad = NULL, 
                   report = TRUE,
                   ref = NULL){
  
  # two things can happen:
  # either Z and beta are lists when different regressions for each off diagonal element
  if(is.list(Z) & is.list(beta)){
    
    # check inputs
    if(length(Z) != length(beta)){
      stop("'Z 'and 'beta' need to have the same length if they are lists.")
    }
    if(!all(sapply(Z, is.matrix))){
      stop("'Z' and 'beta' need to be lists of matrices.")
    }
    # if(!any(sapply(beta, is.vector) | sapply(beta, is.matrix))){
    #   stop("Some entries in 'beta' are neither vectors nor matrices.")
    # }
    beta <- lapply(beta, as.matrix)
    if(any(sapply(Z, ncol) != sapply(beta, nrow))){
      stop("Some entries in 'Z' and 'beta' have a different number of columns.")
    }
    if(any(sapply(Z, nrow) != nrow(Z[[1]]))){
      stop("All matrices in 'Z' need to have the same number of rows.")
    }
    if(any(sapply(beta, ncol) != 1)){
      stop("Matrices in 'beta' can only have one column.")
    }
    
    K <- length(beta)
    nObs <- nrow(Z[[1]])
    Eta <- matrix(NaN, nrow = nObs, ncol = K)
    
    for(i in seq_along(Z)){
      Eta[, i] <- Z[[i]] %*% beta[[i]] # linear predictor matrix
    }
    
  } else if(is.matrix(beta)){
    
    if(is.list(Z)){
      stop("'Z' needs to be a vector or matrix if 'beta' is a matrix.")
    } else if(is.vector(Z)){
      Z <- as.matrix(Z)
    }
    
    K <- nrow(beta)
    p <- ncol(beta) - 1
    
    if(ncol(Z) == p){
      Z = cbind(1, Z) # adding intercept column
    } else if(ncol(Z) != p + 1){
      stop("The dimensions of 'Z' and 'beta' do not match.")
    }
    
    Eta <- Z %*% t(beta) # linear predictor matrix
    
    colnames(beta) <- colnames(Z) # nicer naming
  }
  
  # computing the number of associated states
  N <- as.integer(0.5 + sqrt(0.25 + K), 0)
  
  if(is.null(ref)){
    ref <- 1:N
  } else{
    if(length(ref) != N){
      stop("The reference state vector 'ref' needs to have length equal to the number of states.")
    }
  }
  
  # report quantities for easy use later
  if(report) {
    # Setting rownames: depends on byrow
    names <- outer(paste0("S", 1:N, ">"), paste0("S", 1:N), FUN = paste0) # matrix
    names[cbind(1:N, ref)] <- NA
    names <- na.omit(if (byrow) c(t(names)) else c(names))
    
    # naming beta based on names and colnames of Z
    if(is.list(beta)){
      names(beta) <- names
      for(i in 1:K) rownames(beta[[i]]) <- colnames(Z[[i]])
    } else {
      rownames(beta) <- names
      colnames(beta) <- colnames(Z)
    }
    
    REPORT(beta)
  }
  
  # if ad is not explicitly provided, check if delta is an advector
  if(is.null(ad)){
    
    if(is.list(beta)){
      ad <- any(sapply(beta, inherits, what = "advector"))
    } else{
      ad <- inherits(beta, "advector")
    }
  }
  
  if(!ad) {
    
    Gamma = tpm_g2_cpp(Eta, N, byrow, ref) # C++ version
    
  } else if(ad) {
    "[<-" <- ADoverload("[<-") # overloading assignment operators, currently necessary
    "c" <- ADoverload("c")
    # "diag<-" <- ADoverload("diag<-")
    
    expEta = exp(Eta)
    Gamma = array(1, dim = c(N, N, nrow(expEta)))
    
    ## Loop over entries (stuff over time happens vectorised which speeds up the tape)
    col_ind <- 1
    for(i in seq_len(N)){ # loop over rows
      for(j in seq_len(N)){ # loop over columns
        if(i != ref[j]){ # only if non-diagonal
          if(byrow){
            Gamma[i, j, ] <- expEta[, col_ind]
          } else{
            Gamma[j, i, ] <- expEta[, col_ind]
          }
          # increase col_ind by one
          col_ind = col_ind + 1
        }
      }
    }
    
    # Normalise rows to sum to 1
    for(i in seq_len(N)){
      # transposing is necessary because Gamma[i,,] / colSums(Gamma[i,,]) does not work as expected
      Gamma[i, , ] <- t(t(Gamma[i, , ]) / rowSums(t(Gamma[i, , ])))
    }
  }
  
  # naming
  statenames <- paste0("S", 1:N)
  rownames(Gamma) <- statenames
  colnames(Gamma) <- statenames
  
  Gamma
}


#' Compute the transition probability matrix of a thinned periodically inhomogeneous Markov chain.
#'
#' If the transition probability matrix of an inhomogeneous Markov chain varies only periodically (with period length \eqn{L}), it converges to a so-called periodically stationary distribution. 
#' This happens, because the thinned Markov chain, which has a full cycle as each time step, has homogeneous transition probability matrix
#' \deqn{\Gamma_t = \Gamma^{(t)} \Gamma^{(t+1)} \dots \Gamma^{(t+L-1)}} for all \eqn{t = 1, \dots, L.}
#' This function calculates the matrix above efficiently as a preliminery step to calculating the periodically stationary distribution.
#' 
#' @param Gamma array of transition probability matrices of dimension c(N,N,L).
#' @param t integer index of the time point in the cycle, for which to calculate the thinned transition probility matrix
#'
#' @return thinned transition probabilty matrix of dimension c(N,N)
#' @export
#'
#' @examples
#' # setting parameters for trigonometric link
#' beta = matrix(c(-1, -2, 2, -1, 2, -4), nrow = 2, byrow = TRUE)
#' # calculating periodically varying t.p.m. array (of length 24 here)
#' Gamma = tpm_p(beta = beta)
#' # calculating t.p.m. of thinned Markov chain
#' tpm_thinned(Gamma, 4)
tpm_thinned = function(Gamma, t){
  tpm_thinned_t_cpp(Gamma, t)
}


#' Build all transition probability matrices of a periodically inhomogeneous HMM
#'
#' @description
#' Given a periodically varying variable such as time of day or day of year and the associated cycle length, 
#' this function calculates the transition probability matrices by applying the inverse multinomial logistic link (also known as softmax) to linear predictors of the form
#' \deqn{ 
#'  \eta^{(t)}_{ij} = \beta_0^{(ij)} + \sum_{k=1}^K \bigl( \beta_{1k}^{(ij)} \sin(\frac{2 \pi k t}{L}) + \beta_{2k}^{(ij)} \cos(\frac{2 \pi k t}{L}) \bigr) }
#' for the off-diagonal elements (\eqn{i \neq j}) of the transition probability matrix.
#' This is relevant for modeling e.g. diurnal variation and the flexibility can be increased by adding smaller frequencies (i.e. increasing \eqn{K}).
#' 
#' @details
#' Note that using this function inside the negative log-likelihood function is convenient, but it performs the basis expansion into sine and cosine terms each time it is called. 
#' As these do not change during the optimisation, using \code{\link{tpm_g}} with a pre-calculated (by \code{\link{trigBasisExp}}) design matrix would be more efficient.
#'
#' @family transition probability matrix functions
#'
#' @param tod equidistant sequence of a cyclic variable
#' 
#' For time of day and e.g. half-hourly data, this could be 1, ..., L and L = 48, or 0.5, 1, 1.5, ..., 24 and L = 24.
#' @param L length of one full cycle, on the scale of tod
#' 
#' @param beta matrix of coefficients for the off-diagonal elements of the transition probability matrix
#' 
#' Needs to be of dimension c(N *(N-1), 2*degree+1), where the first column contains the intercepts.
#' @param degree degree of the trigonometric link function
#' 
#' For each additional degree, one sine and one cosine frequency are added.
#' @param Z pre-calculated design matrix (excluding intercept column)
#' 
#' Defaults to \code{NULL} if trigonometric link should be calculated. 
#' From an efficiency perspective, \code{Z} should be pre-calculated within the likelihood function, as the basis expansion should not be redundantly calculated. This can be done by using \code{\link{trigBasisExp}}.
#' @param byrow logical indicating if each transition probability matrix should be filled by row
#'  
#' Defaults to \code{FALSE}, but should be set to \code{TRUE} if one wants to work with a matrix of \code{beta} parameters returned by popular HMM packages like \code{moveHMM}, \code{momentuHMM}, or \code{hmmTMB}.
#' @param ad optional logical, indicating whether automatic differentiation with RTMB should be used. By default, the function determines this itself.
#' @param report logical, indicating whether the coefficient matrix \code{beta} should be reported from the fitted model. Defaults to \code{TRUE}, but only works if \code{ad = TRUE}.
#'
#' @return array of transition probability matrices of dimension c(N,N,length(tod))
#' @export
#'
#' @examples
#' # hourly data 
#' tod = seq(1, 24, by = 1)
#' L = 24
#' beta = matrix(c(-1, 2, -1, -2, 1, -1), nrow = 2, byrow = TRUE)
#' Gamma = tpm_p(tod, L, beta, degree = 1)
#' 
#' # half-hourly data
#' ## integer tod sequence
#' tod = seq(1, 48, by = 1)
#' L = 48
#' beta = matrix(c(-1, 2, -1, -2, 1, -1), nrow = 2, byrow = TRUE)
#' Gamma1 = tpm_p(tod, L, beta, degree = 1)
#' 
#' ## equivalent specification
#' tod = seq(0.5, 24, by = 0.5)
#' L = 24
#' beta = matrix(c(-1, 2, -1, -2, 1, -1), nrow = 2, byrow = TRUE)
#' Gamma2 = tpm_p(tod, L, beta, degree = 1)
#' 
#' all(Gamma1 == Gamma2) # same result
tpm_p = function(tod = 1:24, L=24, beta, degree = 1, Z = NULL, byrow = FALSE, ad = NULL, report = TRUE){
  K = nrow(beta)
  p = ncol(beta) - 1 # number of covariates
  # for N > 1: K = N*(N-1) is bijective with solution
  N = as.integer(0.5 + sqrt(0.25+K), 0)
  
  # check wheter design matrix is provided
  # if not, build trigonometric design matrix
  if(is.null(Z)){
    Z = cbind(1, trigBasisExp(tod, L, degree))
  } else{
    if(ncol(Z) == p){ # intercept column is missing
      Z = cbind(1, Z) # adding intercept column
    } else if(ncol(Z) != p + 1){
      stop("The dimensions of Z and beta do not match.")
    }
  }
  
  # just an interface to tpm_g
  tpm_g(Z, beta, byrow, ad, report)
}



# Functions for continuous-time models ------------------------------------


#' Calculate continuous time transition probabilities
#'
#' A continuous-time Markov chain is described by an infinitesimal generator matrix \eqn{Q}. 
#' When observing data at time points \eqn{t_1, \dots, t_n} the transition probabilites between \eqn{t_i} and \eqn{t_{i+1}} are caluclated as
#' \deqn{\Gamma(\Delta t_i) = \exp(Q \Delta t_i),}
#' where \eqn{\exp()} is the matrix exponential. The mapping \eqn{\Gamma(\Delta t)} is also called the \strong{Markov semigroup}.
#' This function calculates all transition matrices based on a given generator and time differences.
#' 
#' @family transition probability matrix functions
#'
#' @param Q infinitesimal generator matrix of the continuous-time Markov chain of dimension c(N,N)
#' @param timediff time differences between observations of length n-1 when based on n observations
#' @param ad optional logical, indicating whether automatic differentiation with \code{RTMB} should be used. By default, the function determines this itself.
#' @param report logical, indicating whether \code{Q} should be reported from the fitted model. Defaults to \code{TRUE}, but only works if \code{ad = TRUE}.
#'
#' @return array of continuous-time transition matrices of dimension c(N,N,n-1)
#' @export
#' @import RTMB
#'
#' @examples
#' # building a Q matrix for a 3-state cont.-time Markov chain
#' Q = generator(rep(-2, 6))
#'
#' # draw random time differences
#' timediff = rexp(100, 10)
#'
#' # compute all transition matrices
#' Gamma = tpm_cont(Q, timediff)
tpm_cont = function(Q, timediff, ad = NULL, report = TRUE){
  
  # report quantities for easy use later
  if(report) {
    RTMB::REPORT(Q)
  }
  
  # if ad is not explicitly provided, check if delta is an advector
  if(is.null(ad)){
    # check if Q has any of the allowed classes
    if(!any(class(Q) %in% c("advector", "numeric", "matrix", "array"))){
      stop("Q needs to be either a vector, matrix or advector.")
    }
    
    # if Q is advector, run ad version of the function
    ad = inherits(Q, "advector")
  }
  
  if(!ad) {
    
    Qube = semigroup_cpp(Q, timediff) # C++ version
    
  } else if(ad) { # ad version with RTMB
    "[<-" <- ADoverload("[<-") # currently necessary
    
    n = length(timediff)
    N = nrow(Q)
    
    Qube = array(NaN, dim = c(N, N, n))
    
    for(t in 1:n){
      Qube[,,t] = as.matrix(Matrix::expm(Q * timediff[t])) # Matrix::expm for AD
    }
  }
  
  # naming
  N <- nrow(Q)
  statenames <- paste0("S", 1:N)
  rownames(Qube) <- statenames
  colnames(Qube) <- statenames
  
  Qube
}


#' Build the generator matrix of a continuous-time Markov chain
#' 
#' This function builds the \strong{infinitesimal generator matrix} for a \strong{continuous-time Markov chain} from an unconstrained parameter vector.
#' 
#' @family transition probability matrix functions
#' 
#' @param param unconstrained parameter vector of length N*(N-1) where N is the number of states of the Markov chain
#' @param byrow logical indicating if the transition probability matrix should be filled by row
#' @param report logical, indicating whether the generator matrix Q should be reported from the fitted model. Defaults to \code{TRUE}, but only works if when automatic differentiation with \code{RTMB} is used.
#'
#' @return infinitesimal generator matrix of dimension c(N,N)
#' @export
#' @import RTMB
#'
#' @examples
#' # 2 states: 2 free off-diagonal elements
#' generator(rep(-1, 2))
#' # 3 states: 6 free off-diagonal elements
#' generator(rep(-2, 6))
generator = function(param, byrow = FALSE, report = TRUE) {
  
  "[<-" <- ADoverload("[<-") # currently necessary
  
  K = length(param)
  # for N > 1: N*(N-1) is bijective with solution
  N = as.integer(0.5 + sqrt(0.25 + K), 0)
  
  Q = diag(N)
  Q[!Q] = exp(param)
  diag(Q) = 0
  
  if(byrow) Q = t(Q) # transpose if necessary
  
  diag(Q) = -rowSums(Q)
  
  if(report) {
    RTMB::REPORT(Q)
  }
  
  # naming
  statenames <- paste0("S", 1:N)
  rownames(Q) <- statenames
  colnames(Q) <- statenames
  
  Q
}



# Functions for HSMMs -----------------------------------------------------


#' Build the embedded transition probability matrix of an HSMM from unconstrained parameter vector
#'
#' @description
#' Hidden semi-Markov models are defined in terms of state durations and an \strong{embedded} transition probability matrix that contains the conditional transition probabilities given that the \strong{current state is left}. This matrix necessarily has diagonal entries all equal to zero as self-transitions are impossible.
#' 
#' This function builds such an embedded/ conditional transition probability matrix from an unconstrained parameter vector. 
#' For each row of the matrix, the inverse multinomial logistic link is applied.
#' 
#' For a matrix of dimension c(N,N), the number of free off-diagonal elements is N*(N-2), hence also the length of \code{param}.
#' This means, for 2 states, the function needs to be called without any arguments, for 3-states with a vector of length 3, for 4 states with a vector of length 8, etc.
#'
#' Compatible with automatic differentiation by \code{RTMB}
#' 
#' @family transition probability matrix functions
#'
#' @param param unconstrained parameter vector of length N*(N-2) where N is the number of states of the Markov chain
#'
#' If the function is called without \code{param}, it will return the conditional transition probability matrix for a 2-state HSMM, which is fixed with 0 diagonal entries and off-diagonal entries equal to 1.
#'
#' @return embedded/ conditional transition probability matrix of dimension c(N,N)
#' @export
#' @import RTMB
#'
#' @examples
#' # 2 states: no free off-diagonal elements
#' omega = tpm_emb()
#' 
#' # 3 states: 3 free off-diagonal elements
#' param = rep(0, 3)
#' omega = tpm_emb(param)
#' 
#' # 4 states: 8 free off-diagonal elements
#' param = rep(0, 8)
#' omega = tpm_emb(param)
tpm_emb = function(param = NULL){
  "[<-" <- ADoverload("[<-") # overloading assignment operators, currently necessary
  "c" <- ADoverload("c")
  "diag<-" <- ADoverload("diag<-")
  
  if(is.null(param)){
    omega = matrix(c(0,1,1,0), nrow = 2, ncol = 2)
  } else{
    K = length(param)
    # for N > 2: N*(N-2) is bijective with solution
    N = as.integer(1 + sqrt(1 + K), 0)
    
    omega = matrix(0,N,N)
    omega[!diag(N)] = as.vector(t(matrix(c(rep(1,N), exp(param)), N, N-1)))
    omega = t(omega) / apply(omega, 2, sum)
  }
  
  # naming
  N = nrow(omega)
  statenames <- paste0("S", 1:N)
  rownames(omega) <- statenames
  colnames(omega) <- statenames
  
  omega
}

#' Build all embedded transition probability matrices of an inhomogeneous HSMM
#'
#' @description
#' Hidden semi-Markov models are defined in terms of state durations and an \strong{embedded} transition probability matrix that contains the conditional transition probabilities given that the \strong{current state is left}. This matrix necessarily has diagonal entries all equal to zero as self-transitions are impossible.
#' We can allow this matrix to vary with covariates, which is the purpose of this function.
#'
#' It builds all embedded/ conditional transition probability matrices based on a design and parameter matrix.
#' For each row of the matrix, the inverse multinomial logistic link is applied.
#' 
#' For a matrix of dimension c(N,N), the number of free off-diagonal elements is N*(N-2) which determines the number of rows of the parameter matrix.
#'
#' Compatible with automatic differentiation by \code{RTMB}
#' 
#' @family transition probability matrix functions
#' 
#' @param Z covariate design matrix with or without intercept column, i.e. of dimension c(n, p) or c(n, p+1)
#'
#' If \code{Z} has only p columns, an intercept column of ones will be added automatically.
#' @param beta matrix of coefficients for the off-diagonal elements of the embedded transition probability matrix
#'
#' Needs to be of dimension c(N*(N-2), p+1), where the first column contains the intercepts.
#' p can be 0, in which case the model is homogeneous.
#' @param report logical, indicating whether the coefficient matrix beta should be reported from the fitted model. Defaults to \code{TRUE}.
#'
#' @return array of embedded/ conditional transition probability matrices of dimension c(N,N,n)
#' @export
#' @import RTMB
#'
#' @examples
#' ## parameter matrix for 3-state HSMM
#' beta = matrix(c(rep(0, 3), -0.2, 0.2, 0.1), nrow = 3)
#' # no intercept
#' Z = rnorm(100)
#' omega = tpm_emb_g(Z, beta)
#' # intercept
#' Z = cbind(1, Z)
#' omega = tpm_emb_g(Z, beta)
tpm_emb_g = function(Z, beta, report = TRUE){
  "[<-" <- ADoverload("[<-") # overloading assignment operators, currently necessary
  "c" <- ADoverload("c")
  "diag<-" <- ADoverload("diag<-")
  
  beta = as.matrix(beta) # turn beta into matrix if it is a vector
  # that way, we can have homogeneity as a special case with Z = rep(1, ...)
  
  p = ncol(beta) - 1 # number of parameters per state (excluding intercept)
  K = nrow(beta) # number of rows in beta -> N*(N-2) off-diagonal elements
  # for N > 2: N*(N-2) is bijective with solution
  N = as.integer(1 + sqrt(1 + K))
  
  Z = as.matrix(Z)
  
  if(ncol(Z) == p){
    Z = cbind(1, Z) # adding intercept column
  } else if(ncol(Z) != p + 1){
    stop("The dimensions of Z and beta do not match.")
  }
  
  n = nrow(Z)
  
  if(report){
    RTMB::REPORT(beta)
  }
  
  expEta = exp(Z %*% t(beta))
  
  omega = array(NaN, dim = c(N,N,n))
  
  for(t in 1:n){
    O = matrix(0, N, N)
    O[!diag(N)] = as.vector(t(matrix(c(rep(1, N), expEta[t,]), N, N-1)))
    omega[,,t] = t(O) / apply(O, 2, sum)
  }
  
  # naming
  statenames <- paste0("S", 1:N)
  rownames(omega) <- statenames
  colnames(omega) <- statenames
  
  omega
}


#' Builds the transition probability matrix of an HSMM-approximating HMM
#'
#' @description
#' Hidden semi-Markov models (HSMMs) are a flexible extension of HMMs, where the state duration distribution is explicitly modelled.
#' For direct numerical maximum likelhood estimation, HSMMs can be represented as HMMs on an enlarged state space (of size \eqn{M}) and with structured transition probabilities.
#' 
#' This function computes the transition matrix to approximate a given HSMM by an HMM with a larger state space.
#'
#' @param omega embedded transition probability matrix of dimension c(N,N) as computed by \code{\link{tpm_emb}}.
#' @param dm state dwell-time distributions arranged in a list of length(N). Each list element needs to be a vector of length N_i, where N_i is the state aggregate size.
#' @param Fm optional list of length N containing the cumulative distribution functions of the dwell-time distributions.
#' @param sparse logical, indicating whether the output should be a \strong{sparse} matrix. Defaults to \code{TRUE}.
#' @param eps rounding value: If an entry of the transition probabily matrix is smaller, than it is rounded to zero. Usually, this should not be changed.
#'
#' @return extended-state-space transition probability matrix of the approximating HMM
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
tpm_hsmm <- function(omega, dm, 
                     Fm = NULL, sparse = TRUE, eps = 1e-10) {
  "[<-" <- ADoverload("[<-")
  "c" <- ADoverload("c")
  "diag<-" <- ADoverload("diag<-")
  
  Nv = sapply(dm, length)
  N = length(Nv)
  
  # Pre-allocate G matrix
  total_cols = sum(Nv)
  G = matrix(0, total_cols, total_cols)
  
  # compute cdf if not provided
  if(is.null(Fm)) Fm = lapply(1:N, function(i) cumsum(c(0, dm[[i]][-Nv[i]])))
  
  row_start = 1  # Track row start index for G
  for (i in 1:N) {
    Ni = Nv[i]
    # ci = dm[[i]] / (1 - Fm[[i]] + eps)
    ci = max2(dm[[i]], eps) / (1 - Fm[[i]] + eps/2)
    cim = max2(1 - ci, 0)
    
    Gi = matrix(0, Ni, total_cols)  # Pre-allocate the block for Gi
    
    col_start = 1  # Track column start index for Gi
    for (j in 1:N) {
      Nj = Nv[j]
      
      if (i == j) {
        if (Ni == 1) {
          Gi[1, (col_start + Nj - 1)] = cim
        } else {
          diag_block = diag(cim[-Ni], Ni - 1, Ni - 1)
          Gi[1:(Ni - 1), (col_start + 1):(col_start + Nj - 1)] = diag_block
          Gi[Ni, (col_start + Nj - 1)] = cim[Ni]
        }
      } else {
        if (Ni == 1) {
          Gi[1, col_start:(col_start + Nj - 1)] = c(omega[i, j] * (1-cim), rep(0, Nj - 1))
        } else {
          Gi[, col_start] = omega[i, j] * (1-cim)
        }
      }
      
      col_start = col_start + Nj
    }
    
    G[row_start:(row_start + Ni - 1), ] = Gi
    row_start = row_start + Ni
  }
  if(sparse){
    G = methods::as(G, "sparseMatrix")
  }
  G
}

#' Builds all transition probability matrices of an inhomogeneous-HSMM-approximating HMM
#'
#' @description
#' Hidden semi-Markov models (HSMMs) are a flexible extension of HMMs. For direct numerical maximum likelhood estimation, HSMMs can be represented as HMMs on an enlarged state space (of size \eqn{M}) and with structured transition probabilities.
#' 
#' This function computes the transition matrices of a periodically inhomogeneos HSMMs.
#'
#' @param omega embedded transition probability matrix
#'
#' Either a matrix of dimension c(N,N) for homogeneous conditional transition probabilities (as computed by \code{\link{tpm_emb}}), or an array of dimension c(N,N,n) for inhomogeneous conditional transition probabilities (as computed by \code{\link{tpm_emb_g}}).
#' @param dm state dwell-time distributions arranged in a list of length N
#'
#' Each list element needs to be a matrix of dimension c(n, N_i), where each row t is the (approximate) probability mass function of state i at time t.
#' @param eps rounding value: If an entry of the transition probabily matrix is smaller, than it is rounded to zero. Usually, this should not be changed.
#'
#' @return list of dimension length \code{n - max(sapply(dm, ncol))}, containing sparse extended-state-space transition probability matrices for each time point (except the first \code{max(sapply(dm, ncol)) - 1}).
#' @export
#'
#' @examples
#' N = 2
#' # time-varying mean dwell times
#' n = 100
#' z = runif(n)
#' beta = matrix(c(2, 2, 0.1, -0.1), nrow = 2)
#' Lambda = exp(cbind(1, z) %*% t(beta))
#' sizes = c(15, 15) # approximating chain with 30 states
#' # state dwell-time distributions
#' dm = lapply(1:N, function(i) sapply(1:sizes[i]-1, dpois, lambda = Lambda[,i]))
#' 
#' ## homogeneous conditional transition probabilites
#' # diagonal elements are zero, rowsums are one
#' omega = matrix(c(0,1,1,0), nrow = N, byrow = TRUE)
#' 
#' # calculating extended-state-space t.p.m.s
#' Gamma = tpm_ihsmm(omega, dm)
#' 
#' ## inhomogeneous conditional transition probabilites
#' # omega can be an array
#' omega = array(omega, dim = c(N,N,n))
#' 
#' # calculating extended-state-space t.p.m.s
#' Gamma = tpm_ihsmm(omega, dm)
tpm_ihsmm = function(omega, dm, 
                     eps = 1e-10){
  n = nrow(dm[[1]]) # length of timeseries
  N = length(dm) # number of states
  Ni = sapply(dm, ncol) # state aggregate sizes
  bigN = sum(Ni) # overall number of approximating states
  maxNi = max(Ni) # maximum state aggregate size
  
  ## computing all cdfs
  Fm = lapply(1:N, function(i) cbind(0, t(apply(dm[[i]][,-Ni[i]], 1, cumsum))))
  
  ## if omega is just matrix, turn into array
  if(is.matrix(omega)){
    omega = array(omega, dim = c(N, N, n))
  }
  
  # for each timepoint, use tpm_hsmm applied to shifted pmfs and cdfs
  Gamma = vector("list", n - maxNi + 1)
  for(k in 1:(n - maxNi + 1)){
    t = k + maxNi - 1
    dmt = Fmt = vector("list", N)
    
    for(i in 1:N){
      ind = t:(t - Ni[i] + 1)
      
      dmt[[i]] = diag(dm[[i]][ind,])
      Fmt[[i]] = diag(Fm[[i]][ind,])
    }
    # apply tpm_hsmm() to shifted pmfs and pdfs
    Gamma[[k]] = tpm_hsmm(omega[,,t], dmt, Fm = Fmt, sparse = TRUE, eps = eps)
  }
  Gamma
}

#' Builds all transition probability matrices of an periodic-HSMM-approximating HMM
#'
#' @description
#' Hidden semi-Markov models (HSMMs) are a flexible extension of HMMs. For direct numerical maximum likelhood estimation, HSMMs can be represented as HMMs on an enlarged state space (of size \eqn{M}) and with structured transition probabilities.
#' 
#' This function computes the transition matrices of a periodically inhomogeneos HSMMs.
#'
#' @param omega embedded transition probability matrix
#'
#' Either a matrix of dimension c(N,N) for homogeneous conditional transition probabilities (as computed by \code{\link{tpm_emb}}), or an array of dimension c(N,N,L) for inhomogeneous conditional transition probabilities (as computed by \code{\link{tpm_emb_g}}).
#' @param dm state dwell-time distributions arranged in a list of length N
#'
#' Each list element needs to be a matrix of dimension c(L, N_i), where each row t is the (approximate) probability mass function of state i at time t.
#' @param eps rounding value: If an entry of the transition probabily matrix is smaller, than it is rounded to zero. Usually, this should not be changed.
#'
#' @return list of dimension length L, containing sparse extended-state-space transition probability matrices of the approximating HMM for each time point of the cycle.
#' @export
#'
#' @examples
#' N = 2 # number of states
#' L = 24 # cycle length
#' # time-varying mean dwell times
#' Z = trigBasisExp(1:L) # trigonometric basis functions design matrix
#' beta = matrix(c(2, 2, 0.1, -0.1, -0.2, 0.2), nrow = 2)
#' Lambda = exp(cbind(1, Z) %*% t(beta))
#' sizes = c(20, 20) # approximating chain with 40 states
#' # state dwell-time distributions
#' dm = lapply(1:N, function(i) sapply(1:sizes[i]-1, dpois, lambda = Lambda[,i]))
#' 
#' ## homogeneous conditional transition probabilites
#' # diagonal elements are zero, rowsums are one
#' omega = matrix(c(0,1,1,0), nrow = N, byrow = TRUE)
#' 
#' # calculating extended-state-space t.p.m.s
#' Gamma = tpm_phsmm(omega, dm)
#' 
#' ## inhomogeneous conditional transition probabilites
#' # omega can be an array
#' omega = array(omega, dim = c(N,N,L))
#' 
#' # calculating extended-state-space t.p.m.s
#' Gamma = tpm_phsmm(omega, dm)
tpm_phsmm = function(omega, dm, 
                     eps = 1e-10){
  L = nrow(dm[[1]]) # period length
  N = length(dm) # number of states
  Ni = sapply(dm, ncol) # state aggregate sizes
  bigN = sum(Ni) # overall number of approximating states
  
  ## computing all cdfs
  Fm = lapply(1:N, function(i) cbind(0, t(apply(dm[[i]][,-Ni[i]], 1, cumsum))))
  
  ## if omega is just matrix, turn into array
  if(is.matrix(omega)){
    omega = array(omega, dim = c(N, N, L))
  }
  
  # for each timepoint, use tpm_hsmm applied to shifted pmfs and cdfs
  Gamma = vector("list", L)
  
  for(t in 1:L){
    dmt = Fmt = vector("list", N)
    
    for(i in 1:N){
      ind = (t:(t - Ni[i] + 1)) %% L
      ind[ind == 0] = L
      
      dmt[[i]] = diag(dm[[i]][ind,])
      Fmt[[i]] = diag(Fm[[i]][ind,])
    }
    # apply tpm_hsmm() to shifted pmfs and pdfs
    Gamma[[t]] = tpm_hsmm(omega[,,t], dmt, Fm = Fmt, sparse = TRUE, eps = eps)
  }
  Gamma
}

