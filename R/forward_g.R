#' General \href{https://www.taylorfrancis.com/books/mono/10.1201/b20790/hidden-markov-models-time-series-walter-zucchini-iain-macdonald-roland-langrock}{forward algorithm} with time-varying transition probability matrix
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
#' ## generating data from inhomogeneous 2-state HMM
#' mu = c(0, 6)
#' sigma = c(2, 4)
#' beta = matrix(c(-2,-2,0.5,-0.5),nrow=2)
#' delta = c(0.5, 0.5)
#' # simulation
#' n = 2000
#' s = x = rep(NA, n)
#' z = rnorm(n, 0, 2)
#' s[1] = sample(1:2, 1, prob = delta)
#' x[1] = rnorm(1, mu[s[1]], sigma[s[1]])
#' for(t in 2:n){
#'   Gamma = diag(2)
#'   Gamma[!Gamma] = exp(beta[,1]+beta[,2]*z[t])
#'   Gamma = Gamma / rowSums(Gamma)
#'   s[t] = sample(1:2, 1, prob = Gamma[s[t-1],])
#'   x[t] = rnorm(1, mu[s[t]], sigma[s[t]])
#' }
#' 
#' ## negative log likelihood function
#' mllk = function(theta.star, x, z){
#'   # parameter transformations for unconstraint optimization
#'   beta = matrix(theta.star[1:4], 2, 2)
#'   Gamma = tpm_g(Z = z, beta = beta)
#'   delta = c(plogis(theta.star[5]), 1-plogis(theta.star[5]))
#'   mu = theta.star[6:7]
#'   sigma = exp(theta.star[8:9])
#'   # calculate all state-dependent probabilities
#'   allprobs = matrix(1, length(x), 2)
#'   for(j in 1:2){ allprobs[,j] = dnorm(x, mu[j], sigma[j]) }
#'   # return negative for minimization
#'   -forward_g(delta, Gamma, allprobs)
#' }
#' 
#' ## fitting an HMM to the data
#' theta.star = c(-2,-2,1,-1,0,0,5,log(2),log(3))
#' mod = nlm(mllk, theta.star, x = x, z = z)
#'
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
  }
  
  # if ad is not explicitly provided, check if delta is an advector
  if(is.null(ad)){
    # check if delta has any of the allowed classes
    if(!any(class(delta) %in% c("advector", "numeric", "matrix", "array"))){
      stop("delta needs to be either a vector, matrix or advector.")
    }
    
    # if delta is advector, run ad version of the function
    ad = inherits(delta, "advector")
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
