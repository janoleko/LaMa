

# HMMs, ct-HMMs and MMPPs -------------------------------------------------


#' \href{https://www.taylorfrancis.com/books/mono/10.1201/b20790/hidden-markov-models-time-series-walter-zucchini-iain-macdonald-roland-langrock}{Forward algorithm} with homogeneous transition probability matrix
#' 
#' Calculates the log-likelihood of a sequence of observations under a homogeneous hidden Markov model using the \strong{forward algorithm}.
#' 
#' @family forward algorithms
#'
#' @param delta initial or stationary distribution of length N, or matrix of dimension c(k,N) for k independent tracks, if \code{trackID} is provided
#' @param Gamma transition probability matrix of dimension c(N,N), or array of k transition probability matrices of dimension c(N,N,k), if \code{trackID} is provided
#' @param allprobs matrix of state-dependent probabilities/ density values of dimension c(n, N)
#' @param trackID optional vector of length n containing IDs
#' 
#' If provided, the total log-likelihood will be the sum of each track's likelihood contribution.
#' In this case, \code{Gamma} can be a matrix, leading to the same transition probabilities for each track, or an array of dimension c(N,N,k), with one (homogeneous) transition probability matrix for each track.
#' Furthermore, instead of a single vector \code{delta} corresponding to the initial distribution, a \code{delta} matrix of initial distributions, of dimension c(k,N), can be provided, such that each track starts with it's own initial distribution.
#' @param ad optional logical, indicating whether automatic differentiation with \code{RTMB} should be used. By default, the function determines this itself.
#' @param report logical, indicating whether \code{delta}, \code{Gamma} and \code{allprobs} should be reported from the fitted model. Defaults to \code{TRUE}, but only works if \code{ad = TRUE}.
#'
#' @return log-likelihood for given data and parameters
#' @export
#' @import RTMB
#'
#' @examples
#' ## negative log likelihood function
#' nll = function(par, step) {
#'  # parameter transformations for unconstrained optimisation
#'  Gamma = tpm(par[1:2]) # multinomial logit link
#'  delta = stationary(Gamma) # stationary HMM
#'  mu = exp(par[3:4])
#'  sigma = exp(par[5:6])
#'  # calculate all state-dependent probabilities
#'  allprobs = matrix(1, length(step), 2)
#'  ind = which(!is.na(step))
#'  for(j in 1:2) allprobs[ind,j] = dgamma2(step[ind], mu[j], sigma[j])
#'  # simple forward algorithm to calculate log-likelihood
#'  -forward(delta, Gamma, allprobs)
#' }
#' 
#' ## fitting an HMM to the trex data
#' par = c(-2,-2,            # initial tpm params (logit-scale)
#'         log(c(0.3, 2.5)), # initial means for step length (log-transformed)
#'         log(c(0.2, 1.5))) # initial sds for step length (log-transformed)
#' mod = nlm(nll, par, step = trex$step[1:1000])
forward = function(delta, Gamma, allprobs, 
                   trackID = NULL, ad = NULL, report = TRUE){
  
  # report quantities for easy use later
  if(report) {
    RTMB::REPORT(delta)
    RTMB::REPORT(Gamma)
    RTMB::REPORT(allprobs)
    if(!is.null(trackID)){
      RTMB::REPORT(trackID)
    }
    type = "homogeneous"
    REPORT(type)
  }
  
  # if ad is not explicitly provided, check if delta is an advector
  if(is.null(ad)){
    # check if delta has any of the allowed classes
    if(!any(class(delta) %in% c("advector", "numeric", "matrix", "array"))){
      stop("delta needs to be either a vector, matrix or advector.")
    }
    
    # if delta is advector, run ad version of the function
    ad = inherits(delta, "advector") | inherits(Gamma, "advector")
  }
  
  # non-ad version in C++
  if(!ad) {
    
    if(inherits(delta, "advector")){
      stop("It seems you forgot to set ad = TRUE.")
    }
    
    if(is.null(trackID)) {
      l = forward_cpp_h(allprobs, delta, Gamma)
    } else {
      trackInd = calc_trackInd(trackID)
      
      k = length(trackInd) # number of tracks
      
      if(is.vector(delta)){
        delta = matrix(delta, nrow = k, ncol = length(delta), byrow = TRUE)
      } else if(dim(delta)[1] != k){
        stop("Delta needs to be either a vector of length N or a matrix of dimension c(k,N), matching the number tracks.")
      }
      
      if(length(dim(Gamma)) == 3){
        if(dim(Gamma)[3]!= k){
          stop("Gamma needs to be either a matrix of dimension c(N,N) or an array of dimension c(N,N,k), matching the number tracks.")
        }
      } else if(is.matrix(Gamma)){
        Gamma = array(Gamma, dim = c(dim(Gamma), k))
      } else{
        stop("Gamma needs to be either a matrix of dimension c(N,N) or an array of dimension c(N,N,k), matching the number tracks.")
      }
      
      l = forward_cpp_h_tracks(allprobs, delta, Gamma, trackInd)
    }
  } else if(ad) { # ad version
    
    "[<-" <- RTMB::ADoverload("[<-") # overloading assignment operators, currently necessary
    "c" <- ADoverload("c")
    "diag<-" <- ADoverload("diag<-")
    
    # if(report) {
    #   RTMB::REPORT(delta)
    #   RTMB::REPORT(Gamma)
    #   RTMB::REPORT(allprobs)
    # }
    
    if(is.null(trackID)) {
      
      delta = matrix(delta, nrow = 1, ncol = length(delta), byrow = TRUE) # reshape to matrix
      
      # forward algorithm
      foo = delta %*% RTMB::diag(allprobs[1,])
      sumfoo = sum(foo)
      phi = foo / sumfoo
      l = log(sumfoo)
      
      for(t in 2:nrow(allprobs)) {
        foo = phi %*% Gamma %*% RTMB::diag(allprobs[t,])
        sumfoo = sum(foo)
        phi = foo / sumfoo
        l = l + log(sumfoo)
      }
      
    } else if(!is.null(trackID)) {
      
      # RTMB::REPORT(trackID)
      
      uID = unique(trackID) # unique track IDs
      k = length(uID) # number of tracks
      N = ncol(allprobs) # number of states
      
      ## dealing with the initial distribution, either a vector of length N 
      # or a matrix of dimension c(k,N) for k tracks
      # if(is.vector(delta)){
      #   Delta = matrix(delta, nrow = k, ncol = N, byrow = TRUE)
      # } else if(is.matrix(delta)){
      #   if(nrow(delta) == 1){
      #     Delta = matrix(c(delta), nrow = k, ncol = N, byrow = TRUE)
      #   } else if(nrow(delta) == k){
      #     Delta = delta
      #   } else {
      #     stop("Delta needs to be either a vector of length N or a matrix of dimension c(k,N), matching the number tracks.")
      #   }
      # }
      
      delta = as.matrix(delta) # reshape to matrix for easier handling
      
      if(ncol(delta) != N) delta = t(delta) # transpose if necessary
      
      if(nrow(delta) == 1) {
        Delta = matrix(delta, nrow = k, ncol = N, byrow = TRUE)
      } else {
        Delta = delta
      }
      
      ## dealing with Gamma, 
      # either a matrix of dimension c(N,N) or an array of dimension c(N,N,k)
      
      if(length(dim(Gamma)) == 3) {
        if(dim(Gamma)[3] != k) stop("Gamma needs to be either a matrix of dimension c(N,N) or an array of dimension c(N,N,k), matching the number tracks.")
      } else if(is.matrix(Gamma)) {
        Gamma = array(Gamma, dim = c(dim(Gamma), k))
      } else {
        stop("Gamma needs to be either a matrix of dimension c(N,N) or an array of dimension c(N,N,k), matching the number tracks.")
      }
      
      ## forward algorithm
      l = 0 # initialize log-likelihood
      for(i in 1:k) {
        ind = which(trackID == uID[i]) # indices of track i
        
        deltai = RTMB::matrix(Delta[i,], nrow = 1, ncol = N)
        
        foo = deltai %*% RTMB::diag(allprobs[ind[1],])
        sumfoo = sum(foo)
        phi = foo / sumfoo
        l = l + log(sumfoo)
        
        Gamma_i = Gamma[,,i]
        
        for(t in 2:length(ind)) {
          foo = phi %*% Gamma_i %*% RTMB::diag(allprobs[ind[t],])
          sumfoo = sum(foo)
          phi = foo / sumfoo
          l = l + log(sumfoo)
        }
      }
    }
  }
  
  return(l)
}


#' General \href{https://www.taylorfrancis.com/books/mono/10.1201/b20790/hidden-markov-models-time-series-walter-zucchini-iain-macdonald-roland-langrock}{forward algorithm} with time-varying transition probability matrix
#'
#' Calculates the log-likelihood of a sequence of observations under a hidden Markov model with time-varying transition probabilities using the \strong{forward algorithm}.
#' 
#' @family forward algorithms
#'
#' @param delta initial or stationary distribution of length N, or matrix of dimension c(k,N) for k independent tracks, if \code{trackID} is provided
#' @param Gamma array of transition probability matrices of dimension c(N,N,n-1), as in a time series of length n, there are only n-1 transitions. 
#' 
#' If an array of dimension c(N,N,n) for a single track is provided, the first slice will be ignored.
#'  
#' If the elements of \eqn{\Gamma^{(t)}} depend on covariate values at t or covariates t+1 is your choice in the calculation of the array, prior to using this function.
#' When conducting the calculation by using tpm_g(), the choice comes down to including the covariate matrix Z[-1,] oder Z[-n,].
#' 
#' If trackInd is provided, Gamma needs to be an array of dimension c(N,N,n), matching the number of rows of allprobs. For each track, the transition matrix at the beginning will be ignored.
#' If the parameters for Gamma are pooled across tracks or not, depends on your calculation of Gamma. If pooled, you can use tpm_g(Z, beta) to calculate the entire array of transition matrices when Z is of dimension c(n,p). \cr
#' 
#' This function can also be used to fit continuous-time HMMs, where each array entry is the Markov semigroup \eqn{\Gamma(\Delta t) = \exp(Q \Delta t)} and \eqn{Q} is the generator of the continuous-time Markov chain.
#' 
#' @param allprobs matrix of state-dependent probabilities/ density values of dimension c(n, N)
#' @param trackID optional vector of length n containing IDs
#' 
#' If provided, the total log-likelihood will be the sum of each track's likelihood contribution.
#' In this case, \code{Gamma} needs to be an array of dimension c(N,N,n), matching the number of rows of allprobs. For each track, the transition matrix at the beginning of the track will be ignored (as there is no transition between tracks).
#' Furthermore, instead of a single vector \code{delta} corresponding to the initial distribution, a \code{delta} matrix of initial distributions, of dimension c(k,N), can be provided, such that each track starts with it's own initial distribution.
#' @param ad optional logical, indicating whether automatic differentiation with \code{RTMB} should be used. By default, the function determines this itself.
#' @param report logical, indicating whether \code{delta}, \code{Gamma} and \code{allprobs} should be reported from the fitted model. Defaults to \code{TRUE}, but only works if \code{ad = TRUE}.
#'
#' @return log-likelihood for given data and parameters
#' @export
#' @import RTMB
#'
#' @examples
#' ## Simple usage
#' Gamma = array(c(0.9, 0.2, 0.1, 0.8), dim = c(2,2,10))
#' delta = c(0.5, 0.5)
#' allprobs = matrix(0.5, 10, 2)
#' forward_g(delta, Gamma, allprobs)
#' 
#' \donttest{
#' ## Full model fitting example
#' ## negative log likelihood function
#' nll = function(par, step, Z) {
#'  # parameter transformations for unconstrained optimisation
#'  beta = matrix(par[1:6], nrow = 2)
#'  Gamma = tpm_g(Z, beta) # multinomial logit link for each time point
#'  delta = stationary(Gamma[,,1]) # stationary HMM
#'  mu = exp(par[7:8])
#'  sigma = exp(par[9:10])
#'  # calculate all state-dependent probabilities
#'  allprobs = matrix(1, length(step), 2)
#'  ind = which(!is.na(step))
#'  for(j in 1:2) allprobs[ind,j] = dgamma2(step[ind], mu[j], sigma[j])
#'  # simple forward algorithm to calculate log-likelihood
#'  -forward_g(delta, Gamma, allprobs)
#' }
#' 
#' ## fitting an HMM to the trex data
#' par = c(-1.5,-1.5,        # initial tpm intercepts (logit-scale)
#'         rep(0, 4),        # initial tpm slopes
#'         log(c(0.3, 2.5)), # initial means for step length (log-transformed)
#'         log(c(0.2, 1.5))) # initial sds for step length (log-transformed)
#' mod = nlm(nll, par, step = trex$step[1:500], Z = trigBasisExp(trex$tod[1:500]))
#' }
forward_g = function(delta, Gamma, allprobs, 
                     trackID = NULL, ad = NULL, report = TRUE) {
  
  # report quantities for easy use later
  if(report) {
    RTMB::REPORT(delta)
    RTMB::REPORT(Gamma)
    RTMB::REPORT(allprobs)
    if(!is.null(trackID)){
      RTMB::REPORT(trackID)
    }
    type = "inhomogeneous"
    REPORT(type)
  }
  
  # if ad is not explicitly provided, check if delta is an advector
  if(is.null(ad)){
    # check if delta has any of the allowed classes
    if(!any(class(delta) %in% c("advector", "numeric", "matrix", "array"))){
      stop("delta needs to be either a vector, matrix or advector.")
    }
    
    # if delta is advector, run ad version of the function
    ad = inherits(delta, "advector") | inherits(Gamma, "advector")
  }
  
  if(!ad) {
    
    n = nrow(allprobs) # number of observations
    
    if(is.null(trackID)){ # only one track
      
      if(dim(Gamma)[3]==n){
        Gamma = Gamma[,,-1]
      }
      l = forward_cpp_g(allprobs, delta, Gamma)
      
    } else{ # several tracks
      
      trackInd = calc_trackInd(trackID) # creating trackInd for faster C++
      
      k = length(trackInd) # number of tracks
      
      if(dim(Gamma)[3]!=n) stop("Gamma needs to be an array of dimension c(N,N,n), matching the number of rows of allprobs.")
      
      if(is.vector(delta)){
        delta = matrix(delta, nrow = k, ncol = length(delta), byrow = TRUE)
      }
      
      l = forward_cpp_g_tracks(allprobs, delta, Gamma, trackInd)
    }
    
  } else if(ad) {
    
    "[<-" <- RTMB::ADoverload("[<-") # overloading assignment operators, currently necessary
    "c" <- ADoverload("c")
    "diag<-" <- ADoverload("diag<-")
    
    if(report) { # report these quantities by default
      RTMB::REPORT(delta)
      RTMB::REPORT(Gamma)
      RTMB::REPORT(allprobs)
    }
    
    N = ncol(allprobs) # number of states
    n = nrow(allprobs) # number of observations
    
    if(is.null(trackID)) { # only one track
      
      delta = matrix(delta, nrow = 1, ncol = N, byrow = TRUE) # reshape to matrix
      
      Gamma = array(Gamma, dim = dim(Gamma))
      if(dim(Gamma)[3] == n) Gamma = Gamma[,,-1] # deleting first slice
      
      # forward algorithm
      foo = delta * allprobs[1,]
      sumfoo = sum(foo)
      phi = foo / sumfoo
      l = log(sumfoo)
      
      for(t in 2:n) {
        foo = (phi %*% Gamma[,,t-1]) * allprobs[t,]
        sumfoo = sum(foo)
        phi = foo / sumfoo
        l = l + log(sumfoo)
      }
      
    } else if(!is.null(trackID)) {
      
      RTMB::REPORT(trackID)
      
      uID = unique(trackID) # unique track IDs
      k = length(uID) # number of tracks
      
      ## dealing with the initial distribution, either a vector of length N 
      # or a matrix of dimension c(k,N) for k tracks
      
      # if(is.vector(delta)){
      #   Delta = matrix(delta, nrow = k, ncol = N, byrow = TRUE)
      # } else if(is.matrix(delta)){
      #   if(nrow(delta) == 1){
      #     Delta = matrix(c(delta), nrow = k, ncol = N, byrow = TRUE)
      #   } else if(nrow(delta) == k){
      #     Delta = delta
      #   } else {
      #     stop("Delta needs to be either a vector of length N or a matrix of dimension c(k,N), matching the number tracks.")
      #   }
      # }
      
      delta = as.matrix(delta) # reshape to matrix for easier handling
      
      if(ncol(delta) != N) delta = t(delta) # transpose if necessary
      
      if(nrow(delta) == 1) {
        Delta = matrix(delta, nrow = k, ncol = N, byrow = TRUE)
      } else {
        Delta = delta
      }
      
      ## dealing with Gamma, needs to be array of dimension c(N,N,k)
      
      if(dim(Gamma)[3] != n) stop("Gamma needs to be an array of dimension c(N,N,n), matching the number of rows of allprobs.")
      Gamma = array(Gamma, dim = dim(Gamma))
      
      ## forward algorithm
      l = 0 # initialize log-likelihood
      for(i in 1:k) {
        ind = which(trackID == uID[i]) # indices of track i
        
        deltai = matrix(Delta[i,], nrow = 1, ncol = N)
        
        foo = deltai * allprobs[ind[1],]
        sumfoo = sum(foo)
        phi = foo / sumfoo
        l = l + log(sumfoo)
        
        for(t in 2:length(ind)) {
          foo = (phi %*% Gamma[,,ind[t]]) * allprobs[ind[t],]
          sumfoo = sum(foo)
          phi = foo / sumfoo
          l = l + log(sumfoo)
        }
      }
    }
  }
  
  return(l)
}


#' \href{https://www.taylorfrancis.com/books/mono/10.1201/b20790/hidden-markov-models-time-series-walter-zucchini-iain-macdonald-roland-langrock}{Forward algorithm} with for periodically varying transition probability matrices
#'
#' Calculates the log-likelihood of a sequence of observations under a hidden Markov model with periodically varying transition probabilities using the \strong{forward algorithm}.
#' 
#' @family forward algorithms
#' 
#' @details
#' When the transition probability matrix only varies periodically (e.g. as a function of time of day), there are only \eqn{L} unique matrices if \eqn{L} is the period length (e.g. \eqn{L=24} for hourly data and time-of-day variation).
#' Thus, it is much more efficient to only calculate these \eqn{L} matrices and index them by a time variable (e.g. time of day or day of year) instead of calculating such a matrix for each index in the data set (which would be redundant).
#' This function allows for that by only expecting a transition probability matrix for each time point in a period and an integer valued (\eqn{1, \dots, L}) time variable that maps the data index to the according time.
#'
#' @param delta initial or stationary distribution of length N, or matrix of dimension c(k,N) for k independent tracks, if \code{trackID} is provided
#' @param Gamma array of transition probability matrices of dimension c(N,N,L).
#' 
#' Here we use the definition \eqn{\Pr(S_t=j \mid S_{t-1}=i) = \gamma_{ij}^{(t)}} such that the transition probabilities between time point \eqn{t-1} and \eqn{t} are an element of \eqn{\Gamma^{(t)}}.
#' @param allprobs matrix of state-dependent probabilities/ density values of dimension c(n, N)
#' @param tod (Integer valued) variable for cycle indexing in 1, ..., L, mapping the data index to a generalised time of day (length n)
#' 
#' For half-hourly data L = 48. It could, however, also be day of year for daily data and L = 365.
#' @param trackID optional vector of length n containing IDs
#' 
#' If provided, the total log-likelihood will be the sum of each track's likelihood contribution.
#' Instead of a single vector \code{delta} corresponding to the initial distribution, a \code{delta} matrix of initial distributions of dimension c(k,N), can be provided, such that each track starts with it's own initial distribution.
#' @param ad optional logical, indicating whether automatic differentiation with \code{RTMB} should be used. By default, the function determines this itself.
#' @param report logical, indicating whether \code{delta}, \code{Gamma} and \code{allprobs} should be reported from the fitted model. Defaults to \code{TRUE}, but only works if \code{ad = TRUE}.
#'
#'
#' @return log-likelihood for given data and parameters
#' @export
#' @importFrom RTMB REPORT
#'
#' @examples
#' ## negative log likelihood function
#' nll = function(par, step, tod) {
#'  # parameter transformations for unconstrained optimisation
#'  beta = matrix(par[1:6], nrow = 2)
#'  Gamma = tpm_p(1:24, beta = beta) # multinomial logit link for each time point
#'  delta = stationary_p(Gamma, tod[1]) # stationary HMM
#'  mu = exp(par[7:8])
#'  sigma = exp(par[9:10])
#'  # calculate all state-dependent probabilities
#'  allprobs = matrix(1, length(step), 2)
#'  ind = which(!is.na(step))
#'  for(j in 1:2) allprobs[ind,j] = dgamma2(step[ind], mu[j], sigma[j])
#'  # simple forward algorithm to calculate log-likelihood
#'  -forward_p(delta, Gamma, allprobs, tod)
#' }
#' 
#' ## fitting an HMM to the nessi data
#' par = c(-2,-2,            # initial tpm intercepts (logit-scale)
#'         rep(0, 4),        # initial tpm slopes
#'         log(c(0.3, 2.5)), # initial means for step length (log-transformed)
#'         log(c(0.2, 1.5))) # initial sds for step length (log-transformed)
#' mod = nlm(nll, par, step = trex$step[1:500], tod = trex$tod[1:500])
forward_p = function(delta, Gamma, allprobs, tod, trackID = NULL, ad = NULL, report = TRUE){
  utod = unique(tod) # unique time of days
  L = length(utod) # cycle length
  
  if(dim(Gamma)[3] != L){
    stop("Gamma needs to be an array of dimension c(N,N,L), matching the number of unique time points in tod.") 
  } 
  
  Gammanew = Gamma[,,tod]
  
  if(report){
    RTMB::REPORT(tod)
    type = "periodic"
    REPORT(type)
  }
  
  forward_g(delta, Gammanew, allprobs, 
            trackID = trackID, ad = ad, report = report)
}




# HSMMs -------------------------------------------------------------------


#' \href{https://www.taylorfrancis.com/books/mono/10.1201/b20790/hidden-markov-models-time-series-walter-zucchini-iain-macdonald-roland-langrock}{Forward algorithm} for homogeneous hidden semi-Markov models
#' 
#' Calculates the (approximate) log-likelihood of a sequence of observations under a homogeneous hidden semi-Markov model using a modified \strong{forward algorithm}.
#' 
#' @family forward algorithms
#'
#' @details
#' Hidden semi-Markov models (HSMMs) are a flexible extension of HMMs, where the state duration distribution is explicitly modelled by a distribution on the positive integers.
#' For direct numerical maximum likelhood estimation, HSMMs can be represented as HMMs on an enlarged state space (of size \eqn{M}) and with structured transition probabilities.
#' 
#' This function is designed to be used with automatic differentiation based on the \code{R} package \code{RTMB}. It will be very slow without it!
#'
#' @param dm list of length N containing vectors of dwell-time probability mass functions (PMFs) for each state. The vector lengths correspond to the approximating state aggregate sizes, hence there should be little probablity mass not covered by these.
#' @param omega matrix of dimension c(N,N) of conditional transition probabilites, also called embedded transition probability matrix. 
#' 
#' Contains the transition probabilities given that the current state is left. Hence, the diagonal elements need to be zero and the rows need to sum to one. Can be constructed using \code{\link{tpm_emb}}.
#' @param allprobs matrix of state-dependent probabilities/ density values of dimension c(n, N) which will automatically be converted to the appropriate dimension.
#' @param trackID optional vector of length n containing IDs
#' 
#' If provided, the total log-likelihood will be the sum of each track's likelihood contribution.
#' In this case, \code{dm} can be a nested list, where the top layer contains k \code{dm} lists as described above. \code{omega} can then also be an array of dimension c(N,N,k) with one conditional transition probability matrix for each track.
#' Furthermore, instead of a single vector \code{delta} corresponding to the initial distribution, a \code{delta} matrix of initial distributions, of dimension c(k,N), can be provided, such that each track starts with it's own initial distribution.
#' @param delta optional vector of initial state probabilities of length N
#' 
#' By default, the stationary distribution is computed (which is typically recommended).
#' @param eps small value to avoid numerical issues in the approximating transition matrix construction. Usually, this should not be changed.
#' @param report logical, indicating whether initial distribution, approximating transition probability matrix and \code{allprobs} matrix should be reported from the fitted model. Defaults to \code{TRUE}.
#'
#' @return log-likelihood for given data and parameters
#' @export
#' @import RTMB
#'
#' @examples
#' # currently no examples
forward_hsmm <- function(dm, omega, allprobs,
                         trackID = NULL, delta = NULL, eps = 1e-10, report = TRUE){
  # overloading assignment operators, currently necessary
  "[<-" <- ADoverload("[<-")
  "c" <- ADoverload("c")
  "diag<-" <- ADoverload("diag<-")
  ################################
  
  agsizes = vapply(dm, length, 1L)
  
  N = ncol(allprobs) # number of HSMM states
  M = sum(agsizes) # total number of states of the approximating HMM
  
  stationary = is.null(delta) # if delta is not provided, stationary distribution needs to be computed
  
  if(is.null(trackID)) {
    ## compute approximating tpm
    Gamma_sparse = tpm_hsmm(omega, dm, eps = eps)
    # Gamma = Matrix::Matrix(Gamma_sparse, sparse = FALSE) # dense for reporting
    
    ## if stationary, compute initial stationary distribution
    if(stationary){
      delta = stationary_sparse(Gamma_sparse)
      delta_sparse = methods::as(t(delta), "sparseMatrix")
    } else{ # if delta is provided, "stuff out" with zeros
      cols_to_fill = c(1, cumsum(agsizes[-N])+1)
      delta_sparse = Matrix::sparseMatrix(i = rep(1, N), j = cols_to_fill, 
                                          x = delta, dims = c(1, M))
    }
    
    ## report quantities for state decoding
    if(report){
      RTMB::REPORT(delta_sparse)
      RTMB::REPORT(Gamma_sparse)
      RTMB::REPORT(allprobs)
    }
    
    # forward algorithm
    # foo = delta_sparse %*% Matrix::Diagonal(x = rep(allprobs[1,], times = agsizes))
    foo = delta_sparse %*% diag(rep(allprobs[1,], times = agsizes))
    sumfoo = sum(foo)
    phi = foo / sumfoo
    l = log(sumfoo)
    
    for(t in 2:nrow(allprobs)) {
      # foo = phi %*% Gamma_sparse %*% Matrix::Diagonal(x = rep(allprobs[t,], times = agsizes))
      foo = phi %*% Gamma_sparse %*% diag(rep(allprobs[t,], times = agsizes))
      sumfoo = sum(foo)
      phi = foo / sumfoo
      l = l + log(sumfoo)
    }
    
  } else if(!is.null(trackID)) {
    
    RTMB::REPORT(trackID) # report trackID for viterbi etc.
    
    uID = unique(trackID) # unique track IDs
    k = length(uID) # number of tracks
    
    if(length(dm) == k){ # dm for each track
      if(is.matrix(omega)){ # different dms but same omegas
        omega = array(omega, dim = c(dim(omega), k))
      } else if(length(dim(omega)) == 3){ # different dms and different omegas
        if(dim(omega)[3] != k) stop("omega needs to be either a matrix of dimension c(N,N) or an array of dimension c(N,N,k), matching the number of tracks.")
      }
      
      Gamma_sparse = lapply(1:k, function(k) tpm_hsmm(omega[,,k], dm[[k]], eps = eps)) # build k Gammas based on omega and dm
      
      if(stationary == TRUE){ # delta not provided but computed: each track starts in stationary distribution
        Delta = t(sapply(1:k, function(k) stationary_sparse(Gamma_sparse[[k]]))) # build k deltas based on Gammas
        Delta_sparse = methods::as(Delta, "sparseMatrix")
      } else{
        
        # if delta is provided
        delta = as.matrix(delta) # reshape to matrix for easier handling
        
        if(ncol(delta) != N) delta = t(delta) # transpose if necessary
        
        if(nrow(delta) == 1) {
          Delta = matrix(delta, nrow = k, ncol = N, byrow = TRUE)
        } else {
          Delta = delta
        }
        
        # "stuff out" delta with zeros
        cols_to_fill = c(1, cumsum(agsizes[-N])+1)
        
        row_indices <- rep(1:k, N)
        col_indices <- rep(cols_to_fill, each = k)
        
        values = as.vector(Delta)
        
        Delta_sparse <- Matrix::sparseMatrix(i = row_indices,   # Row indices
                                             j = col_indices,   # Column indices
                                             x = values,        # Non-zero values (from dense matrix)
                                             dims = c(k, M))    # Final sparse matrix dimensions
      }
    } else if(length(dm) == length(agsizes)){ # only one dm
      if(is.matrix(omega)){ # if omega is matrix, just build one Gamma
        Gamma_sparse = tpm_hsmm(omega, dm, eps = eps)
        
        ## if stationary, compute initial stationary distribution
        if(stationary){
          delta = stationary_sparse(Gamma_sparse)
          Delta = matrix(delta, nrow = k, ncol = length(delta), byrow = TRUE) # repeat rows for each track
          Delta_sparse = Matrix::Matrix(Delta, sparse = TRUE) # make sparse
        } else{ # if delta is provided, "stuff out" with zeros
          cols_to_fill = c(1, cumsum(agsizes[-N])+1)
          
          if(is.matrix(delta)){
            values = as.vector(delta)
          } else{
            values = rep(delta, each = k)
          }
          Delta_sparse = Matrix::sparseMatrix(i = rep(1:k, N), j = rep(cols_to_fill, each=k), 
                                              x = values, dims = c(k, M))
        }
        Gamma_sparse = list(Gamma_sparse)
        Gamma_sparse = rep(Gamma_sparse, k)
      } else if(length(dim(omega)) == 3){ # different omegas but same dm
        if(dim(omega)[3] != k) stop("omega needs to be either a matrix of dimension c(N,N) or an array of dimension c(N,N,k), matching the number of tracks.")
        
        Gamma_sparse = lapply(1:k, function(k) tpm_hsmm(omega[,,k], dm, eps = eps)) # build k Gammas with fixed dm but varying omega
        
        if(stationary == TRUE){ # delta not provided but computed: each track starts in stationary distribution
          Delta = t(sapply(1:k, function(k) stationary_sparse(Gamma_sparse[[k]]))) # build k deltas based on Gammas
          Delta_sparse = methods::as(Delta, "sparseMatrix")
        } else{
          # if delta is provided
          delta = as.matrix(delta) # reshape to matrix for easier handling
          
          if(ncol(delta) != N) delta = t(delta) # transpose if necessary
          
          if(nrow(delta) == 1) {
            Delta = matrix(delta, nrow = k, ncol = N, byrow = TRUE)
          } else {
            Delta = delta
          }
          
          # "stuff out" delta with zeros
          cols_to_fill = c(1, cumsum(agsizes[-N])+1)
          
          row_indices <- rep(1:k, N)
          col_indices <- rep(cols_to_fill, each = k)
          
          values = as.vector(Delta)
          
          Delta_sparse <- Matrix::sparseMatrix(i = row_indices,   # Row indices
                                               j = col_indices,   # Column indices
                                               x = values,        # Non-zero values (from dense matrix)
                                               dims = c(k, M))    # Final sparse matrix dimensions
        }
      }
    }
    
    ## report quantities for state decoding
    if(report){
      RTMB::REPORT(Delta_sparse)
      RTMB::REPORT(Gamma_sparse)
      RTMB::REPORT(allprobs)
    }
    
    ## forward algorithm
    l = 0 # initialize log-likelihood
    for(i in 1:k) {
      ind = which(trackID == uID[i]) # indices of track i
      
      delta_i = Delta_sparse[i, , drop = FALSE]
      Gamma_i = Gamma_sparse[[i]]
      
      foo = delta_i %*% Matrix::Diagonal(x = rep(allprobs[ind[1],], times = agsizes))
      sumfoo = sum(foo)
      phi = foo / sumfoo
      l_i = log(sumfoo)
      
      for(t in 2:length(ind)) {
        foo = phi %*% Gamma_i %*% Matrix::Diagonal(x = rep(allprobs[ind[t],], times = agsizes))
        sumfoo = sum(foo)
        phi = foo / sumfoo
        l_i = l_i + log(sumfoo)
      }
      
      l = l + l_i
    }
  }
  
  l
}


#' \href{https://www.taylorfrancis.com/books/mono/10.1201/b20790/hidden-markov-models-time-series-walter-zucchini-iain-macdonald-roland-langrock}{Forward algorithm} for hidden semi-Markov models with inhomogeneous state durations and/ or conditional transition probabilities
#'
#' Calculates the (approximate) log-likelihood of a sequence of observations under an inhomogeneous hidden semi-Markov model using a modified \strong{forward algorithm}.
#' 
#' @family forward algorithms
#' 
#' @details
#' Hidden semi-Markov models (HSMMs) are a flexible extension of HMMs, where the state duration distribution is explicitly modelled by a distribution on the positive integers. This function can be used to fit HSMMs where the state-duration distribution and/ or the conditional transition probabilities vary with covariates.
#' For direct numerical maximum likelhood estimation, HSMMs can be represented as HMMs on an enlarged state space (of size \eqn{M}) and with structured transition probabilities.
#' 
#' This function is designed to be used with automatic differentiation based on the \code{R} package \code{RTMB}. It will be very slow without it!
#'
#' @param dm list of length N containing matrices (or vectors) of dwell-time probability mass functions (PMFs) for each state.
#' 
#' If the dwell-time PMFs are constant, the vectors are the PMF of the dwell-time distribution fixed in time. The vector lengths correspond to the approximating state aggregate sizes, hence there should be little probablity mass not covered by these.
#' 
#' If the dwell-time PMFs are inhomogeneous, the matrices need to have n rows, where n is the number of observations. The number of columns again correponds to the size of the approximating state aggregates.
#' 
#' In the latter case, the first \code{max(sapply(dm, ncol)) - 1} observations will not be used because the first approximating transition probability matrix needs to be computed based on the first \code{max(sapply(dm, ncol))} covariate values (represented by \code{dm}).
#' @param omega matrix of dimension c(N,N) or array of dimension c(N,N,n) of conditional transition probabilites, also called embedded transition probability matrix.
#' 
#' It contains the transition probabilities given the current state is left. Hence, the diagonal elements need to be zero and the rows need to sum to one. Such a matrix can be constructed using \code{\link{tpm_emb}} and an array using \code{\link{tpm_emb_g}}.
#' @param allprobs matrix of state-dependent probabilities/ density values of dimension c(n, N)
#' @param trackID trackID optional vector of length n containing IDs
#' 
#' If provided, the total log-likelihood will be the sum of each track's likelihood contribution.
#' Instead of a single vector \code{delta} corresponding to the initial distribution, a \code{delta} matrix of initial distributions, of dimension c(k,N), can be provided, such that each track starts with it's own initial distribution.
#' @param delta optional vector of initial state probabilities of length N
#' 
#' By default, instead of this, the stationary distribution is computed corresponding to the first approximating transition probability matrix of each track is computed. Contrary to the homogeneous case, this is not theoretically motivated but just for convenience.
#' @param startInd optional integer index at which the forward algorithm starts. 
#' 
#' When approximating inhomogeneous HSMMs by inhomogeneous HMMs, the first transition probability matrix that can be constructed is at time \code{max(sapply(dm, ncol))} (as it depends on the previous covariate values).
#' Hence, when not provided, \code{startInd} is chosen to be \code{max(sapply(dm, ncol))}. Fixing \code{startInd} at a value \strong{larger} than max(aggregate sizes) is useful when models with different aggregate sizes are fitted to the same data and are supposed to be compared. In that case it is important that all models use the same number of observations.
#' @param eps small value to avoid numerical issues in the approximating transition matrix construction. Usually, this should not be changed.
#' @param report logical, indicating whether initial distribution, approximating transition probability matrix and \code{allprobs} matrix should be reported from the fitted model. Defaults to \code{TRUE}.
#'
#' @return log-likelihood for given data and parameters
#' @export
#' @import RTMB
#'
#' @examples
#' # currently no examples
forward_ihsmm <- function(dm, omega, allprobs,
                          trackID = NULL, delta = NULL, startInd = NULL,
                          eps = 1e-10, report = TRUE){
  
  # overloading assignment operators, currently necessary
  "[<-" <- ADoverload("[<-")
  "c" <- ADoverload("c")
  "diag<-" <- ADoverload("diag<-")
  ################################
  
  N = ncol(allprobs) # number of HSMM states
  
  # obtain aggregate sizes based on dm
  if(is.matrix(dm[[1]])){
    agsizes = sapply(dm, ncol)
  } else if(is.vector(dm[[1]])){
    agsizes = sapply(dm, length)
  }
  
  M = sum(agsizes) # total number of states of the approximating HMM
  n = nrow(allprobs) # number of observations
  maxag = max(agsizes) # maximum number of states in the approximating HMM
  # this is important, because the forward algo can only start at maxag because the first 1:(maxag-1) data points are needed for the first tpm
  
  # check if startInd is provided
  # if not, take the maximum aggregate size as startInd
  if(is.null(startInd)){
    startInd = maxag
  } else {
    if(startInd < maxag) stop("startInd must be at least the maximum aggregate size")
  }
  
  stationary = is.null(delta) # if delta is not provided, stationary distribution needs to be computed
  
  if(is.null(trackID)) {
    
    ## first case: dm varies -> truely inhomogeneous HSMM
    if(is.matrix(dm[[1]])){
      if(nrow(dm[[1]]) != n) stop("dm needs to be a list of length N, either containing vectors or matrices of dwell-time PMFs.")
      
      ## compute approximating tpm
      Gamma_sparse = tpm_ihsmm(omega, dm, eps = eps)
      
      ## if stationary, compute initial stationary distribution
      if(stationary){
        delta = stationary_sparse(Gamma_sparse[[startInd - maxag + 1]])
        delta_sparse = methods::as(t(delta), "sparseMatrix")
      } else{ # if delta is provided, "stuff out" with zeros
        cols_to_fill = c(1, cumsum(agsizes[-N])+1)
        delta_sparse = Matrix::sparseMatrix(i = rep(1, N), j = cols_to_fill, 
                                            x = delta, dims = c(1, M))
      }
      
      ## report quantities for state decoding
      if(report) {
        RTMB::REPORT(delta_sparse)
        RTMB::REPORT(Gamma_sparse)
        RTMB::REPORT(allprobs)
      }
      
      # forward algorithm
      # foo = delta_sparse %*% Matrix::Diagonal(x = rep(allprobs[1,], times = agsizes))
      foo = delta_sparse %*% diag(rep(allprobs[startInd,], times = agsizes))
      sumfoo = sum(foo)
      phi = foo / sumfoo
      l = log(sumfoo)
      
      for(t in (startInd + 1):nrow(allprobs)) {
        # foo = phi %*% Gamma_sparse %*% Matrix::Diagonal(x = rep(allprobs[t,], times = agsizes))
        foo = phi %*% Gamma_sparse[[t - maxag]] %*% diag(rep(allprobs[t,], times = agsizes))
        sumfoo = sum(foo)
        phi = foo / sumfoo
        l = l + log(sumfoo)
      }
      
    } else if(is.vector(dm[[1]])){
      if(is.matrix(omega)){ # if omega is matrix, just build one Gamma
        Gamma_sparse = tpm_hsmm(omega, dm, eps = eps)
        
        ## if stationary, compute initial stationary distribution
        if(stationary){
          delta = stationary_sparse(Gamma_sparse)
          delta = matrix(delta, nrow = 1, ncol = length(delta), byrow = TRUE)
          delta_sparse = Matrix::Matrix(delta, sparse = TRUE) # make sparse
        } else{ # if delta is provided, "stuff out" with zeros
          cols_to_fill = c(1, cumsum(agsizes[-N]) + 1)
          delta_sparse = Matrix::sparseMatrix(i = rep(1:k, N), j = rep(cols_to_fill, each=k), 
                                              x = delta, dims = c(k, M))
        }
        Gamma_sparse = list(Gamma_sparse)
        Gamma_sparse = rep(Gamma_sparse, n)
        
      } else if(length(dim(omega)) == 3){ # different omegas but same dm
        if(dim(omega)[3] != n) stop("omega needs to be either a matrix of dimension c(N,N) or an array of dimension c(N,N,n), matching the number of observations")
        
        Gamma_sparse = lapply(1:n, function(k) tpm_hsmm(omega[,,k], dm, eps = eps)) # build k Gammas with fixed dm but varying omega
        
        if(stationary == TRUE){ # delta not provided but computed: each track starts in stationary distribution
          delta = stationary_sparse(Gamma_sparse[[1]]) # build k deltas based on Gammas
          delta_sparse = methods::as(t(delta), "sparseMatrix")
        } else{
          cols_to_fill = c(1, cumsum(agsizes[-N])+1)
          delta_sparse = Matrix::sparseMatrix(i = rep(1, N), j = cols_to_fill, 
                                              x = delta, dims = c(1, M))
        }
      }
      
      ## forward algorithm
      foo = delta_sparse %*% diag(rep(allprobs[1,], times = agsizes))
      sumfoo = sum(foo)
      phi = foo / sumfoo
      l = log(sumfoo)
      
      for(t in 2:nrow(allprobs)) {
        foo = phi %*% Gamma_sparse[[t-1]] %*% diag(rep(allprobs[t,], times = agsizes))
        sumfoo = sum(foo)
        phi = foo / sumfoo
        l = l + log(sumfoo)
      }
    }
    
  } else if(!is.null(trackID)) {
    
    RTMB::REPORT(trackID) # report trackID for viterbi etc.
    
    uID = unique(trackID) # unique track IDs
    k = length(uID) # number of tracks
    
    ## first case: dm varies -> truely inhomogeneous HSMM
    if(is.matrix(dm[[1]])){ # dm varies
      if(nrow(dm[[1]]) != n) stop("dm needs to be a list of length N, either containing vectors or matrices of dwell-time PMFs.")
      if(is.matrix(omega)){ # fixed omega
        omega = array(omega, dim = c(dim(omega), n))
      } else if(length(dim(omega)) == 3){ # different dms and different omegas
        if(dim(omega)[3] != n) stop("omega needs to be either a matrix of dimension c(N,N) or an array of dimension c(N,N,n), matching the number of observations")
      }
      
      # nested list of sparse approximating tpms
      Gamma_sparse = lapply(1:k, function(i){
        ind = which(trackID == uID[i])
        this_dm = lapply(1:N, function(j) dm[[j]][ind,])
        tpm_ihsmm(omega[,,ind], this_dm, eps = eps)
      })
      
      if(stationary == TRUE){ # delta not provided but computed: each track starts in stationary distribution
        trackInd = calc_trackInd(trackID)
        Delta = t(sapply(1:k, function(i) stationary_sparse(Gamma_sparse[[i]][[1]]))) # build k deltas based on Gammas
        Delta_sparse = methods::as(Delta, "sparseMatrix")
      } else{
        
        # if delta is provided
        delta = as.matrix(delta) # reshape to matrix for easier handling
        
        if(ncol(delta) != N) delta = t(delta) # transpose if necessary
        
        if(nrow(delta) == 1) {
          Delta = matrix(delta, nrow = k, ncol = N, byrow = TRUE)
        } else {
          Delta = delta
        }
        
        # "stuff out" delta with zeros
        cols_to_fill = c(1, cumsum(agsizes[-N])+1)
        
        row_indices <- rep(1:k, N)
        col_indices <- rep(cols_to_fill, each = k)
        
        values = as.vector(Delta)
        
        Delta_sparse <- Matrix::sparseMatrix(i = row_indices,   # Row indices
                                             j = col_indices,   # Column indices
                                             x = values,        # Non-zero values (from dense matrix)
                                             dims = c(k, M))    # Final sparse matrix dimensions
      }
      
      ## forward algorithm
      l = 0 # initialize log-likelihood
      for(i in 1:k) {
        ind = which(trackID == uID[i]) # indices of track i
        
        delta_i = Delta_sparse[i, , drop = FALSE]
        Gamma_i = Gamma_sparse[[i]]
        
        # foo = delta_sparse %*% Matrix::Diagonal(x = rep(allprobs[1,], times = agsizes))
        foo = delta_i %*% diag(rep(allprobs[ind[startInd],], times = agsizes))
        sumfoo = sum(foo)
        phi = foo / sumfoo
        l_i = log(sumfoo)
        
        for(t in (startInd + 1):length(ind)) {
          foo = phi %*% Gamma_i[[t-maxag]] %*% diag(rep(allprobs[ind[t],], times = agsizes))
          sumfoo = sum(foo)
          phi = foo / sumfoo
          l_i = l_i + log(sumfoo)
        }
        
        l = l + l_i
      }
      
      ## second case: dm does not vary, only omega varies -> easier
    } else if(is.vector(dm[[1]])){ # only one dm
      if(is.matrix(omega)){ # if omega is matrix, just build one Gamma
        Gamma_sparse = tpm_hsmm(omega, dm, eps = eps)
        
        ## if stationary, compute initial stationary distribution
        if(stationary){
          delta = stationary_sparse(Gamma_sparse)
          Delta = matrix(delta, nrow = k, ncol = length(delta), byrow = TRUE) # repeat rows for each track
          Delta_sparse = Matrix::Matrix(Delta, sparse = TRUE) # make sparse
        } else{ # if delta is provided, "stuff out" with zeros
          cols_to_fill = c(1, cumsum(agsizes[-N]) + 1)
          
          if(is.matrix(delta)){
            values = as.vector(delta)
          } else{
            values = rep(delta, each = k)
          }
          Delta_sparse = Matrix::sparseMatrix(i = rep(1:k, N), j = rep(cols_to_fill, each=k), 
                                              x = values, dims = c(k, M))
        }
        Gamma_sparse = list(Gamma_sparse)
        Gamma_sparse = rep(Gamma_sparse, n)
        
      } else if(length(dim(omega)) == 3){ # different omegas but same dm
        if(dim(omega)[3] != n) stop("omega needs to be either a matrix of dimension c(N,N) or an array of dimension c(N,N,n), matching the number of observations")
        
        Gamma_sparse = lapply(1:n, function(k) tpm_hsmm(omega[,,k], dm, eps = eps)) # build k Gammas with fixed dm but varying omega
        
        if(stationary == TRUE){ # delta not provided but computed: each track starts in stationary distribution
          trackInd = calc_trackInd(trackID)
          Delta = t(sapply(1:k, function(k) stationary_sparse(Gamma_sparse[[trackInd[k]]]))) # build k deltas based on Gammas
          Delta_sparse = methods::as(Delta, "sparseMatrix")
        } else{
          # if delta is provided
          delta = as.matrix(delta) # reshape to matrix for easier handling
          
          if(ncol(delta) != N) delta = t(delta) # transpose if necessary
          
          if(nrow(delta) == 1) {
            Delta = matrix(delta, nrow = k, ncol = N, byrow = TRUE)
          } else {
            Delta = delta
          }
          
          # "stuff out" delta with zeros
          cols_to_fill = c(1, cumsum(agsizes[-N])+1)
          
          row_indices <- rep(1:k, N)
          col_indices <- rep(cols_to_fill, each = k)
          
          values = as.vector(Delta)
          
          Delta_sparse <- Matrix::sparseMatrix(i = row_indices,   # Row indices
                                               j = col_indices,   # Column indices
                                               x = values,        # Non-zero values (from dense matrix)
                                               dims = c(k, M))    # Final sparse matrix dimensions
        }
      }
      
      ## forward algorithm
      l = 0 # initialize log-likelihood
      for(i in 1:k) {
        ind = which(trackID == uID[i]) # indices of track i
        
        delta_i = Delta_sparse[i, , drop = FALSE]
        Gamma_i = Gamma_sparse[ind]
        
        foo = delta_i %*% diag(rep(allprobs[ind[1],], times = agsizes))
        sumfoo = sum(foo)
        phi = foo / sumfoo
        l_i = log(sumfoo)
        
        for(t in 2:length(ind)) {
          foo = phi %*% Gamma_i[[t-1]] %*% diag(rep(allprobs[ind[t],], times = agsizes))
          sumfoo = sum(foo)
          phi = foo / sumfoo
          l_i = l_i + log(sumfoo)
        }
        
        l = l + l_i
      }
    }
  }
  
  l
}


#' \href{https://www.taylorfrancis.com/books/mono/10.1201/b20790/hidden-markov-models-time-series-walter-zucchini-iain-macdonald-roland-langrock}{Forward algorithm} for hidden semi-Markov models with periodically inhomogeneous state durations and/ or conditional transition probabilities
#'
#' Calculates the (approximate) log-likelihood of a sequence of observations under a periodically inhomogeneous hidden semi-Markov model using a modified \strong{forward algorithm}.
#' 
#' @family forward algorithms
#' 
#' @description
#' Hidden semi-Markov models (HSMMs) are a flexible extension of HMMs, where the state duration distribution is explicitly modelled by a distribution on the positive integers. This function can be used to fit HSMMs where the state-duration distribution and/ or the conditional transition probabilities vary with covariates.
#' For direct numerical maximum likelhood estimation, HSMMs can be represented as HMMs on an enlarged state space (of size \eqn{M}) and with structured transition probabilities.
#' 
#' This function can be used to fit HSMMs where the state-duration distribution and/ or the conditional transition probabilities vary periodically.
#' In the special case of periodic variation (as compared to arbitrary covariate influence), this version is to be preferred over \code{\link{forward_ihsmm}} because it computes the \strong{correct periodically stationary distribution} and no observations are lost for the approximation.
#' 
#' This function is designed to be used with automatic differentiation based on the \code{R} package \code{RTMB}. It will be very slow without it!
#'
#' @param dm list of length N containing matrices (or vectors) of dwell-time probability mass functions (PMFs) for each state.
#'
#' If the dwell-time PMFs are constant, the vectors are the PMF of the dwell-time distribution fixed in time. The vector lengths correspond to the approximating state aggregate sizes, hence there should be little probablity mass not covered by these.
#' 
#' If the dwell-time PMFs are inhomogeneous, the matrices need to have L rows, where L is the cycle length. The number of columns again correpond to the size of the approximating state aggregates.
#' @param omega matrix of dimension c(N,N) or array of dimension c(N,N,L) of conditional transition probabilites, also called embedded transition probability matrix
#' 
#' It contains the transition probabilities given the current state is left. Hence, the diagonal elements need to be zero and the rows need to sum to one. Such a matrix can be constructed using \code{\link{tpm_emb}} and an array using \code{\link{tpm_emb_g}}.
#' @param allprobs matrix of state-dependent probabilities/ density values of dimension c(n, N)
#' @param tod (Integer valued) variable for cycle indexing in 1, ..., L, mapping the data index to a generalised time of day (length n).
#' For half-hourly data L = 48. It could, however, also be day of year for daily data and L = 365.
#' @param trackID optional vector of length n containing IDs
#' 
#' If provided, the total log-likelihood will be the sum of each track's likelihood contribution.
#' Instead of a single vector \code{delta} corresponding to the initial distribution, a \code{delta} matrix of initial distributions, of dimension c(k,N), can be provided, such that each track starts with it's own initial distribution.
#' @param delta Optional vector of initial state probabilities of length N. By default, instead of this, the stationary distribution is computed corresponding to the first approximating t.p.m. of each track is computed. Contrary to the homogeneous case, this is not theoretically motivated but just for convenience.
#' @param eps small value to avoid numerical issues in the approximating transition matrix construction. Usually, this should not be changed.
#' @param report logical, indicating whether initial distribution, approximating transition probability matrix and \code{allprobs} matrix should be reported from the fitted model. Defaults to \code{TRUE}.
#'
#' @return log-likelihood for given data and parameters
#' @export
#' @import RTMB
#'
#' @examples
#' # currently no examples
forward_phsmm <- function(dm, omega, allprobs, tod,
                          trackID = NULL, delta = NULL, eps = 1e-10, report = TRUE){
  
  # overloading assignment operators, currently necessary
  "[<-" <- ADoverload("[<-")
  "c" <- ADoverload("c")
  "diag<-" <- ADoverload("diag<-")
  ################################
  
  N = ncol(allprobs) # number of HSMM states
  
  # obtain aggregate sizes based on dm
  if(is.matrix(dm[[1]])){
    agsizes = sapply(dm, ncol)
  } else if(is.vector(dm[[1]])){
    agsizes = sapply(dm, length)
  }
  
  M = sum(agsizes) # total number of states of the approximating HMM
  n = nrow(allprobs) # number of observations
  L = length(unique(tod)) # cycle length -> number of unique tpms
  # this is important, because the forward algo can only start at maxag because the first 1:(maxag-1) data points are needed for the first tpm
  
  stationary = is.null(delta) # if delta is not provided, stationary distribution needs to be computed
  
  if(is.null(trackID)) {
    
    ## first case: dm varies -> truely inhomogeneous HSMM
    if(is.matrix(dm[[1]])){
      if(nrow(dm[[1]]) != L) stop("dm needs to be a list of length N, either containing vectors or matrices of dwell-time PMFs.")
      
      ## compute approximating tpm
      Gamma_sparse = tpm_phsmm(omega, dm, eps = eps)
      
      ## if stationary, compute initial stationary distribution
      if(stationary){
        delta = stationary_p_sparse(Gamma_sparse, t = tod[1])
        delta_sparse = methods::as(t(delta), "sparseMatrix")
      } else{ # if delta is provided, "stuff out" with zeros
        cols_to_fill = c(1, cumsum(agsizes[-N])+1)
        delta_sparse = Matrix::sparseMatrix(i = rep(1, N), j = cols_to_fill, 
                                            x = delta, dims = c(1, M))
      }
      
      ## report quantities for state decoding
      if(report) {
        RTMB::REPORT(delta_sparse)
        RTMB::REPORT(Gamma_sparse)
        RTMB::REPORT(allprobs)
      }
      
      # forward algorithm
      # foo = delta_sparse %*% Matrix::Diagonal(x = rep(allprobs[1,], times = agsizes))
      foo = delta_sparse %*% diag(rep(allprobs[1,], times = agsizes))
      sumfoo = sum(foo)
      phi = foo / sumfoo
      l = log(sumfoo)
      
      for(t in 2:nrow(allprobs)) {
        # foo = phi %*% Gamma_sparse %*% Matrix::Diagonal(x = rep(allprobs[t,], times = agsizes))
        foo = phi %*% Gamma_sparse[[tod[t-1]]] %*% diag(rep(allprobs[t,], times = agsizes))
        sumfoo = sum(foo)
        phi = foo / sumfoo
        l = l + log(sumfoo)
      }
      
    } else if(is.vector(dm[[1]])){
      if(is.matrix(omega)){ # if omega is matrix, just build one Gamma
        Gamma_sparse = tpm_hsmm(omega, dm, eps = eps)
        
        ## if stationary, compute initial stationary distribution
        if(stationary){
          delta = stationary_sparse(Gamma_sparse)
          delta = matrix(delta, nrow = 1, ncol = length(delta), byrow = TRUE)
          delta_sparse = Matrix::Matrix(delta, sparse = TRUE) # make sparse
        } else{ # if delta is provided, "stuff out" with zeros
          cols_to_fill = c(1, cumsum(agsizes[-N]) + 1)
          delta_sparse = Matrix::sparseMatrix(i = rep(1:k, N), j = rep(cols_to_fill, each=k), 
                                              x = delta, dims = c(k, M))
        }
        Gamma_sparse = list(Gamma_sparse)
        Gamma_sparse = rep(Gamma_sparse, L)
        
      } else if(length(dim(omega)) == 3){ # different omegas but same dm
        if(dim(omega)[3] != L) stop("omega needs to be either a matrix of dimension c(N,N) or an array of dimension c(N,N,L), matching the cycle length.")
        
        Gamma_sparse = lapply(1:L, function(k) tpm_hsmm(omega[,,k], dm, eps = eps)) # build k Gammas with fixed dm but varying omega
        
        if(stationary == TRUE){ # delta not provided but computed: each track starts in stationary distribution
          delta = stationary_p_sparse(Gamma_sparse, t = tod[1])
          delta_sparse = methods::as(t(delta), "sparseMatrix")
        } else{
          cols_to_fill = c(1, cumsum(agsizes[-N])+1)
          delta_sparse = Matrix::sparseMatrix(i = rep(1, N), j = cols_to_fill, 
                                              x = delta, dims = c(1, M))
        }
      }
      
      ## forward algorithm
      foo = delta_sparse %*% diag(rep(allprobs[1,], times = agsizes))
      sumfoo = sum(foo)
      phi = foo / sumfoo
      l = log(sumfoo)
      
      for(t in 2:nrow(allprobs)) {
        foo = phi %*% Gamma_sparse[[tod[t-1]]] %*% diag(rep(allprobs[t,], times = agsizes))
        sumfoo = sum(foo)
        phi = foo / sumfoo
        l = l + log(sumfoo)
      }
    }
    
  } else if(!is.null(trackID)) {
    
    RTMB::REPORT(trackID) # report trackID for viterbi etc.
    
    uID = unique(trackID) # unique track IDs
    k = length(uID) # number of tracks
    
    ## first case: dm varies -> truely inhomogeneous HSMM
    if(is.matrix(dm[[1]])){ # dm varies
      if(nrow(dm[[1]]) != L) stop("dm needs to be a list of length N, either containing vectors or matrices of dwell-time PMFs.")
      if(length(dim(omega)) == 3){ # different dms and different omegas
        if(dim(omega)[3] != L) stop("omega needs to be either a matrix of dimension c(N,N) or an array of dimension c(N,N,L), matching the cycle length.")
      }
      
      Gamma_sparse = tpm_phsmm(omega, dm, eps = eps) # still only L unique tpms
      
      if(stationary == TRUE){ # delta not provided but computed: each track starts in stationary distribution
        trackInd = calc_trackInd(trackID)
        Delta = stationary_p_sparse(Gamma_sparse) # calculate all periodically stationaries -> not very expensive because of recursion
        Delta = Delta[tod[trackInd], ] # corresponding initial periodically stationary for each track
        Delta_sparse = methods::as(Delta, "sparseMatrix") 
      } else{
        
        # if delta is provided
        delta = as.matrix(delta) # reshape to matrix for easier handling
        
        if(ncol(delta) != N) delta = t(delta) # transpose if necessary
        
        if(nrow(delta) == 1) {
          Delta = matrix(delta, nrow = k, ncol = N, byrow = TRUE)
        } else {
          Delta = delta
        }
        
        # "stuff out" delta with zeros
        cols_to_fill = c(1, cumsum(agsizes[-N])+1)
        
        row_indices <- rep(1:k, N)
        col_indices <- rep(cols_to_fill, each = k)
        
        values = as.vector(Delta)
        
        Delta_sparse <- Matrix::sparseMatrix(i = row_indices,   # Row indices
                                             j = col_indices,   # Column indices
                                             x = values,        # Non-zero values (from dense matrix)
                                             dims = c(k, M))    # Final sparse matrix dimensions
      }
      
      ## forward algorithm
      l = 0 # initialize log-likelihood
      for(i in 1:k) {
        ind = which(trackID == uID[i]) # indices of track i
        
        delta_i = Delta_sparse[i, , drop = FALSE]
        thistod = tod[ind]
        
        # foo = delta_sparse %*% Matrix::Diagonal(x = rep(allprobs[1,], times = agsizes))
        foo = delta_i %*% diag(rep(allprobs[ind[1],], times = agsizes))
        sumfoo = sum(foo)
        phi = foo / sumfoo
        l_i = log(sumfoo)
        
        for(t in 2:length(ind)) {
          foo = phi %*% Gamma_sparse[[thistod[t-1]]] %*% diag(rep(allprobs[ind[t],], times = agsizes))
          sumfoo = sum(foo)
          phi = foo / sumfoo
          l_i = l_i + log(sumfoo)
        }
        
        l = l + l_i
      }
      
      ## second case: dm does not vary, only omega varies -> easier
    } else if(is.vector(dm[[1]])){ # only one dm
      if(is.matrix(omega)){ # if omega is matrix, just build one Gamma
        Gamma_sparse = tpm_hsmm(omega, dm, eps = eps)
        
        ## if stationary, compute initial stationary distribution
        if(stationary){
          delta = stationary_sparse(Gamma_sparse)
          Delta = matrix(delta, nrow = k, ncol = length(delta), byrow = TRUE) # repeat rows for each track
          Delta_sparse = Matrix::Matrix(Delta, sparse = TRUE) # make sparse
        } else{ # if delta is provided, "stuff out" with zeros
          cols_to_fill = c(1, cumsum(agsizes[-N]) + 1)
          
          if(is.matrix(delta)){
            values = as.vector(delta)
          } else{
            values = rep(delta, each = k)
          }
          Delta_sparse = Matrix::sparseMatrix(i = rep(1:k, N), j = rep(cols_to_fill, each=k), 
                                              x = values, dims = c(k, M))
        }
        Gamma_sparse = list(Gamma_sparse)
        Gamma_sparse = rep(Gamma_sparse, L)
        
      } else if(length(dim(omega)) == 3){ # different omegas but same dm
        if(dim(omega)[3] != L) stop("omega needs to be either a matrix of dimension c(N,N) or an array of dimension c(N,N,L), matching the cycle length.")
        
        Gamma_sparse = lapply(1:L, function(k) tpm_hsmm(omega[,,k], dm, eps = eps)) # build k Gammas with fixed dm but varying omega
        
        if(stationary == TRUE){ # delta not provided but computed: each track starts in stationary distribution
          trackInd = calc_trackInd(trackID)
          Delta = stationary_p_sparse(Gamma_sparse) # calculate all periodically stationaries -> not very expensive because of recursion
          Delta = Delta[tod[trackInd], ] # corresponding initial periodically stationary for each track
          Delta_sparse = methods::as(Delta, "sparseMatrix") 
        } else{
          # if delta is provided
          delta = as.matrix(delta) # reshape to matrix for easier handling
          
          if(ncol(delta) != N) delta = t(delta) # transpose if necessary
          
          if(nrow(delta) == 1) {
            Delta = matrix(delta, nrow = k, ncol = N, byrow = TRUE)
          } else {
            Delta = delta
          }
          
          # "stuff out" delta with zeros
          cols_to_fill = c(1, cumsum(agsizes[-N])+1)
          
          row_indices <- rep(1:k, N)
          col_indices <- rep(cols_to_fill, each = k)
          
          values = as.vector(Delta)
          
          Delta_sparse <- Matrix::sparseMatrix(i = row_indices,   # Row indices
                                               j = col_indices,   # Column indices
                                               x = values,        # Non-zero values (from dense matrix)
                                               dims = c(k, M))    # Final sparse matrix dimensions
        }
      }
      
      ## forward algorithm
      l = 0 # initialize log-likelihood
      for(i in 1:k) {
        ind = which(trackID == uID[i]) # indices of track i
        
        delta_i = Delta_sparse[i, , drop = FALSE]
        thistod = tod[ind]
        
        # foo = delta_sparse %*% Matrix::Diagonal(x = rep(allprobs[1,], times = agsizes))
        foo = delta_i %*% diag(rep(allprobs[ind[1],], times = agsizes))
        sumfoo = sum(foo)
        phi = foo / sumfoo
        l_i = log(sumfoo)
        
        for(t in 2:length(ind)) {
          foo = phi %*% Gamma_sparse[[thistod[t-1]]] %*% diag(rep(allprobs[ind[t],], times = agsizes))
          sumfoo = sum(foo)
          phi = foo / sumfoo
          l_i = l_i + log(sumfoo)
        }
        
        l = l + l_i
      }
    }
  }
  
  l
}
