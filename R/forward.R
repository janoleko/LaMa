#' \href{https://www.taylorfrancis.com/books/mono/10.1201/b20790/hidden-markov-models-time-series-walter-zucchini-iain-macdonald-roland-langrock}{Forward algorithm} with homogeneous transition probability matrix
#'
#' @param delta Initial or stationary distribution of length N, or matrix of dimension c(k,N) for k independent tracks, if trackInd is provided
#' @param Gamma Transition probability matrix of dimension c(N,N), or array of k transition probability matrices of dimension c(N,N,k), if trackInd is provided.
#' @param allprobs Matrix of state-dependent probabilities/ density values of dimension c(n, N)
#' @param trackID Optional vector of length n containing IDs. If provided, the total log-likelihood will be the sum of each track's likelihood contribution.
#' In this case, Gamma can be a matrix, leading to the same transition probabilities for each track, or an array of dimension c(N,N,k), with one (homogeneous) transition probability matrix for each track.
#' Furthermore, instead of a single vector delta corresponding to the initial distribution, a delta matrix of initial distributions, of dimension c(k,N), can be provided, such that each track starts with it's own initial distribution.
#' @param ad Logical, indicating whether automatic differentiation with RTMB should be used. Defaults to FALSE.
#' @param report Logical, indicating whether delta, Gamma and allprobs should be reported from the fitted model. Defaults to TRUE, but only works if ad = TRUE.
#'
#' @return Log-likelihood for given data and parameters
#' @export
#'
#' @examples
#' ## generating data from homogeneous 2-state HMM
#' mu = c(0, 6)
#' sigma = c(2, 4)
#' Gamma = matrix(c(0.5, 0.05, 0.15, 0.85), nrow = 2, byrow = TRUE)
#' delta = c(0.5, 0.5)
#' # simulation
#' s = x = rep(NA, 500)
#' s[1] = sample(1:2, 1, prob = delta)
#' x[1] = rnorm(1, mu[s[1]], sigma[s[1]])
#' for(t in 2:500){
#'   s[t] = sample(1:2, 1, prob = Gamma[s[t-1],])
#'   x[t] = rnorm(1, mu[s[t]], sigma[s[t]])
#' }
#' 
#' ## negative log likelihood function
#' mllk = function(theta.star, x){
#'   # parameter transformations for unconstraint optimization
#'   Gamma = tpm(theta.star[1:2])
#'   delta = stationary(Gamma) # stationary HMM
#'   mu = theta.star[3:4]
#'   sigma = exp(theta.star[5:6])
#'   # calculate all state-dependent probabilities
#'   allprobs = matrix(1, length(x), 2)
#'   for(j in 1:2){ allprobs[,j] = dnorm(x, mu[j], sigma[j]) }
#'   # return negative for minimization
#'   -forward(delta, Gamma, allprobs)
#' }
#' 
#' ## fitting an HMM to the data
#' theta.star = c(-2,-2,0,5,log(2),log(3))
#' mod = stats::nlm(mllk, theta.star, x = x)
#'
forward = function(delta, Gamma, allprobs, 
                   trackID = NULL, ad = FALSE, report = TRUE){
  
  if(!ad) {
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
  } else if(ad) {
    
    "[<-" <- RTMB::ADoverload("[<-")
    
    if(report) {
      RTMB::REPORT(delta)
      RTMB::REPORT(Gamma)
      RTMB::REPORT(allprobs)
    }
    
    if(is.null(trackID)) {
      
      delta = matrix(delta, nrow = 1, ncol = length(delta), byrow = TRUE) # reshape to matrix
      
      # forward algorithm
      foo = delta * allprobs[1,]
      sumfoo = sum(foo)
      phi = foo / sumfoo
      l = log(sumfoo)
      
      for(t in 2:nrow(allprobs)) {
        foo = (phi %*% Gamma) * allprobs[t,]
        sumfoo = sum(foo)
        phi = foo / sumfoo
        l = l + log(sumfoo)
      }
      
    } else if(!is.null(trackID)) {
      
      RTMB::REPORT(trackID)
      
      uID = unique(trackID) # unique track IDs
      k = length(uID) # number of tracks
      N = ncol(allprobs) # number of states
      
      ## dealing with the initial distribution, either a vector of length N 
      # or a matrix of dimension c(k,N) for k tracks
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
        
        deltai = matrix(Delta[i,], nrow = 1, ncol = N)
        
        foo = deltai * allprobs[ind[1],]
        sumfoo = sum(foo)
        phi = foo / sumfoo
        l = l + log(sumfoo)
        
        Gamma_i = Gamma[,,i]
        
        for(t in 2:length(ind)) {
          foo = (phi %*% Gamma_i) * allprobs[ind[t],]
          sumfoo = sum(foo)
          phi = foo / sumfoo
          l = l + log(sumfoo)
        }
      }
    }
  }
  
  return(l)
}
