#' Calculate conditional local state probabilities for inhomogeneous HMMs
#' 
#' Computes \cr \cr
#' \eqn{\Pr(S_t = j \mid X_1, ..., X_T)} \cr \cr
#'
#' @param delta Initial distribution of length N, or matrix of dimension c(k,N) for k independent tracks, if \code{trackID} is provided
#' @param Gamma Array of transition probability matrices of dimension c(N,N,n-1), as in a time series of length n, there are only n-1 transitions. 
#' If you provide an array of dimension c(N,N,n) for a single track, the first slice will be ignored. \cr
#' If you provide \code{trackID}, \code{Gamma} needs to be an array of dimension c(N,N,n), where n is the number of rows in \code{allprobs}. Then for each track the first transition matrix will be ignored.
#' @param allprobs Matrix of state-dependent probabilities/ density values of dimension c(n, N)
#' @param trackID Optional vector of k track IDs, if multiple tracks need to be decoded separately
#'
#' @return Matrix of conditional state probabilities of dimension c(n,N)
#' @export
#'
#' @examples
#' Gamma = tpm_g(runif(99), matrix(c(-1,-1,1,-2), nrow = 2, byrow = TRUE))
#' delta = c(0.5, 0.5)
#' allprobs = matrix(runif(200), nrow = 100, ncol = 2)
#' 
#' probs = stateprobs_g(delta, Gamma, allprobs)
stateprobs_g = function(delta, Gamma, allprobs, trackID = NULL){
  n = nrow(allprobs)
  N = ncol(allprobs)
  
  # If ID is provided, several tracks need to be decoded separately
  if(!is.null(trackID)){
    uID = unique(trackID) # unique IDs
    k = length(uID) # number of tracks
    
    if(is.vector(delta)){ # single initial distribution
      delta = matrix(delta, nrow = k, ncol = length(delta), byrow = TRUE)
    } else if(is.matrix(delta)){
      if(dim(delta)[1] != k){
        if(dim(delta)[1] == 1){
          delta = matrix(c(delta), nrow = k, ncol = length(delta), byrow = TRUE)
        } else{
          stop("Delta needs to be either a vector of length N or a matrix of dimension c(k,N), matching the number tracks.")
        }
      }
    }
    
    if(dim(delta)[2] != N) stop("Initial distribution(s) do not match the number of states implied by allprobs.")
    
    # initialising state vector
    stateprobs = matrix(NA, nrow = n, ncol = N)
    
    # loop over individual tracks
    for(i in 1:length(uID)){
      id_i = which(trackID == uID[i])
      allprobs_i = allprobs[id_i, ]
      
      if(length(id_i) == 1) stop("All tracks must be at least of length 2.")
      
      if(nrow(allprobs_i) == 2){
        # viterbi algorithm for track of length 2 only
        Gamma_i = as.matrix(Gamma[,,id_i[-1]])
        
        lalpha = matrix(NA, 2, N)
        lbeta = matrix(NA, 2, N)
        
        # computing forward variables on log scale
        foo = delta[i,] * allprobs_i[1,]
        l = log(sum(foo))
        foo = foo / sum(foo)
        lalpha[1,] = log(foo) + l
        
        foo = (foo %*% Gamma_i) * allprobs_i[2,]
        l = l + log(sum(foo))
        foo = foo / sum(foo)
        lalpha[2,] = log(foo) + l
        
        # computing backward variables on log scale
        foo = rep(1, N)
        lbeta[2,] = rep(0, N)
        foo = Gamma_i %*% diag(allprobs_i[2, ]) %*% foo
        lbeta[1,] = log(foo)
      
        c = max(lalpha[2, ])
        llk = c + log(sum(exp(lalpha[2, ] - c)))
        stateprobs[id_i, ] = exp(log(alpha) + lbeta - llk)
        
      } else{
        # regurlar local decoding
        lalpha = logalpha_cpp(allprobs_i, delta[i,], Gamma[,,id_i[-1]])
        lbeta = logbeta_cpp(allprobs_i, Gamma[,,id_i[-1]])
        
        c = max(lalpha[nrow(lalpha), ])
        llk = c + log(sum(exp(lalpha[nrow(lalpha), ] - c)))
        
        stateprobs[id_i, ] = exp(lalpha + lbeta - llk)
      }
    }
    
  } else{
    if(!is.vector(delta)){
      stop("If trackID is not provided, delta needs to be a vector of length N.")
    }
    
    if(dim(Gamma)[3] == n){
      warning("Igoring the first slice of Gamma, as there are only n-1 transitions in a time series of length n.")
      # not using the first slice of Gamma, if n slices are provided
      Gamma = Gamma[,,-1]
    }
    
    lalpha = logalpha_cpp(allprobs, delta, Gamma)
    lbeta = logbeta_cpp(allprobs, Gamma)
    
    c = max(lalpha[length(lalpha),])
    llk = c + log(sum(exp(lalpha[length(lalpha),] - c)))
    
    stateprobs = exp(lalpha + lbeta - llk)
  }
  
  # rowSums should already be one, but just to be safe
  stateprobs / rowSums(statprobs)
}
