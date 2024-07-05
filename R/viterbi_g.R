#' Viterbi algorithm for decoding states of inhomogeneous HMMs
#'
#' @param delta Initial distribution of length N, or matrix of dimension c(k,N) for k independent tracks, if ID is provided
#' @param Gamma Array of transition probability matrices of dimension c(N,N,n-1), as in a time series of length n, there are only n-1 transitions. 
#' If you provide an array of dimension c(N,N,n), the first slice will be ignored. \cr
#' If you provide an ID vector, Gamma needs to be an array of dimension c(N,N,n), where n is the number of rows in allprobs. Then for each track the first transition matrix will be ignored.
#' @param allprobs Matrix of state-dependent probabilities/ density values of dimension c(n, N)
#' @param ID Optional vector of k track IDs, if multiple tracks need to be decoded separately
#'
#' @return Vector of decoded states of length n
#' @export
#'
#' @examples
#' delta = c(0.5, 0.5)
#' Gamma = array(dim = c(2,2,99))
#' for(t in 1:99){
#'   gammas = rbeta(2, shape1 = 0.4, shape2 = 1)
#'   Gamma[,,t] = matrix(c(1-gammas[1], gammas[1], 
#'                       gammas[2], 1-gammas[2]), nrow = 2, byrow = TRUE)
#' }
#' allprobs = matrix(runif(200), nrow = 100, ncol = 2)
#' states = viterbi_g(delta, Gamma, allprobs)
viterbi_g = function(delta, Gamma, allprobs, ID = NULL){
  n = nrow(allprobs)
  N = ncol(allprobs)
  
  # If ID is provided, several tracks need to be decoded separately
  if(!is.null(ID)){
    uID = unique(ID)
    k = length(uID) # number of tracks
    
    if(is.vector(delta)){
      delta = matrix(delta, nrow = k, ncol = length(delta), byrow = TRUE)
    } else if(dim(delta)[1] != k){
      stop("Delta needs to be either a vector of length N or a matrix of dimension c(k,N), matching the number tracks.")
    }
    
    # initialising state vector
    states = rep(NA, n)
    
    # loop over individual tracks
    for(i in 1:length(uID)){
      id_i = which(ID == uID[i])
      
      allprobs_i = allprobs[id_i, ]
      
      if(length(id_i) == 1) stop("All tracks must be at least of length 2.")
      
      if(dim(allprobs_i)[1] == 2){
        # viterbi algorithm for track of length 2 only
        Gamma_i = as.matrix(Gamma[,,id_i[-1]])
        
        xi = matrix(0, nrow = 2, ncol = N) 
        foo = delta[i,] * allprobs_i[1, ]
        xi[1, ] = foo / sum(foo)
        
        foo = apply(xi[1, ] * Gamma_i, 2, max) * allprobs[2, ]
        xi[2, ] = foo / sum(foo)
        
        iv = numeric(2)
        iv[2] = which.max(xi[2, ]) 
        iv[1] = which.max(Gamma_i[, iv[2]] * xi[1, ])
        
        states[id_i] = iv
        
      } else{
        # regurlar viterbi algorithm for this track
        states[id_i] = viterbi_g_cpp(allprobs_i, delta[i,], Gamma[,,id_i[-1]])
      }
      
    }
    
  } else{
    if(!is.vector(delta)){
      stop("If no ID is provided, delta needs to be a vector of length N.")
    }
    
    if(dim(Gamma)[3]==n){
      warning("Igoring the first slice of Gamma, as there are only n-1 transitions in a time series of length n.")
      # not using the first slice of Gamma, if n slices are provided
      Gamma = Gamma[,,-1]
    }
    
    states = viterbi_g_cpp(allprobs, delta, Gamma)
  }
  
  return(as.integer(states))
}
  

