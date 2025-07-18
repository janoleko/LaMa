

# Global decoding ---------------------------------------------------------


#' Viterbi algorithm for state decoding in homogeneous HMMs
#' 
#' The Viterbi algorithm allows one to decode the most probable state sequence of an HMM.
#' 
#' @family decoding functions
#'
#' @param delta initial distribution of length N, or matrix of dimension c(k,N) for k independent tracks, if \code{trackID} is provided
#' @param Gamma transition probability matrix of dimension c(N,N) or array of transition probability matrices of dimension c(N,N,k) if \code{trackID} is provided
#' @param allprobs matrix of state-dependent probabilities/ density values of dimension c(n, N)
#' @param trackID optional vector of k track IDs, if multiple tracks need to be decoded separately
#' @param mod optional model object containing initial distribution \code{delta}, transition probability matrix \code{Gamma}, matrix of state-dependent probabilities \code{allprobs}, and potentially a \code{trackID} variable
#' 
#' If you are using automatic differentiation either with \code{RTMB::MakeADFun} or \code{\link{qreml}} and include \code{\link{forward}} in your likelihood function, the objects needed for state decoding are automatically reported after model fitting.
#' Hence, you can pass the model object obtained from running \code{RTMB::report()} or from \code{\link{qreml}} directly to this function.
#'
#' @return vector of decoded states of length n
#' @export
#'
#' @examples
#' delta = c(0.5, 0.5)
#' Gamma = matrix(c(0.9, 0.1, 0.2, 0.8), nrow = 2, byrow = TRUE)
#' allprobs = matrix(runif(200), nrow = 100, ncol = 2)
#' states = viterbi(delta, Gamma, allprobs)
viterbi = function(delta, Gamma, allprobs, trackID = NULL,
                   mod = NULL){
  
  # check if a model with delta, Gamma and allprobs is provided
  if(!is.null(mod)){
    if(is.null(mod$delta)){
      stop("Model object contains no initial distribution.")
    }
    if(is.null(mod$Gamma)){
      stop("Model object contains no transition matrix.")
    }
    if(is.null(mod$allprobs)){
      stop("Model object contains no state-dependent probabilities.")
    }
    
    # if suitable model object is provided, overwrite inputs with model object
    delta = mod$delta
    Gamma = mod$Gamma
    allprobs = mod$allprobs
    trackID = mod$trackID
  }
  
  n = nrow(allprobs)
  N = ncol(allprobs)
  
  # inflating Gamma to use viterbi_g
  if(is.null(trackID)){
    if(length(dim(Gamma)) > 2){
      stop("If no 'trackID' is provided, 'Gamma' needs to be a matrix of dimension c(N,N). \nDo you need to use 'viterbi_g()'?")
    }
    Gammanew = array(Gamma, dim = c(N, N, n-1))
  } else{
    uID = unique(trackID)
    k = length(uID) # number of tracks
    
    if(is.matrix(Gamma)){
      Gammanew = array(Gamma, dim = c(N, N, n))
    } else if(is.array(Gamma)){
      if(dim(Gamma)[3] != k) stop("Number of distinct transition matrices does not match the number of tracks.")
      
      ## construct integer trackID
      integerID = match(trackID, uID)
      
      Gammanew = Gamma[,,integerID]
    }
  }
  
  viterbi_g(delta, Gammanew, allprobs, trackID)
}


#' Viterbi algorithm for state decoding in inhomogeneous HMMs
#' 
#' The Viterbi algorithm allows one to decode the most probable state sequence of an HMM.
#' 
#' @family decoding functions
#' 
#' @param delta initial distribution of length N, or matrix of dimension c(k,N) for k independent tracks, if \code{trackID} is provided
#' @param Gamma array of transition probability matrices of dimension c(N,N,n-1), as in a time series of length n, there are only n-1 transitions
#' 
#' If an array of dimension c(N,N,n) is provided for a single track, the first slice will be ignored.
#' 
#' If \code{trackID} is provided, \code{Gamma} needs to be an array of dimension c(N,N,n), where n is the number of rows in \code{allprobs}. Then for each track the first transition matrix will be ignored.
#' @param allprobs matrix of state-dependent probabilities/ density values of dimension c(n, N)
#' @param trackID optional vector of k track IDs, if multiple tracks need to be decoded separately
#' @param mod optional model object containing initial distribution \code{delta}, transition probability matrix \code{Gamma}, matrix of state-dependent probabilities \code{allprobs}, and potentially a \code{trackID} variable
#' 
#' If you are using automatic differentiation either with \code{RTMB::MakeADFun} or \code{\link{qreml}} and include \code{\link{forward_g}} in your likelihood function, the objects needed for state decoding are automatically reported after model fitting.
#' Hence, you can pass the model object obtained from running \code{RTMB::report()} or from \code{\link{qreml}} directly to this function.
#'
#' @return vector of decoded states of length n
#' @export
#'
#' @examples
#' delta = c(0.5, 0.5)
#' Gamma = tpm_g(runif(10), matrix(c(-2,-2,1,-1), nrow = 2))
#' allprobs = matrix(runif(20), nrow = 10, ncol = 2)
#' states = viterbi_g(delta, Gamma[,,-1], allprobs)
viterbi_g = function(delta, Gamma, allprobs, trackID = NULL,
                     mod = NULL){
  
  # check if a model with delta, Gamma and allprobs is provided
  if(!is.null(mod)){
    if(is.null(mod$delta)){
      stop("Model object contains no initial distribution.")
    }
    if(is.null(mod$Gamma)){
      stop("Model object contains no transition matrix.")
    }
    if(is.null(mod$allprobs)){
      stop("Model object contains no state-dependent probabilities.")
    }
    
    # if suitable model object is provided, overwrite inputs with model object
    delta = mod$delta
    # Gamma = mod$Gamma
    Gamma = mod$Gamma[,, -1] # ignoring first slice
    allprobs = mod$allprobs
    trackID = mod$trackID
  }
  
  n = nrow(allprobs)
  N = ncol(allprobs)
  
  # If ID is provided, several tracks need to be decoded separately
  if(!is.null(trackID)){
    uID = unique(trackID)
    k = length(uID) # number of tracks
    
    if(is.vector(delta)){
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
    
    # initialising state vector
    states = rep(NA, n)
    
    # loop over individual tracks
    for(i in seq_len(length(uID))){
      id_i = which(trackID == uID[i])
      
      allprobs_i = allprobs[id_i, ]
      
      if(length(id_i) == 1) stop("All tracks must be at least of length 2.")
      
      if(dim(allprobs_i)[1] == 2){
        # viterbi algorithm for track of length 2 only
        Gamma_i = as.matrix(Gamma[,,id_i[-1]])
        
        xi = matrix(0, nrow = 2, ncol = N) 
        foo = delta[i,] * allprobs_i[1, ]
        xi[1, ] = foo / sum(foo)
        
        foo = apply(xi[1, ] * Gamma_i, 2, max) * allprobs_i[2, ]
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


#' Viterbi algorithm for state decoding in periodically inhomogeneous HMMs
#' 
#' The Viterbi algorithm allows one to decode the most probable state sequence of an HMM.
#'
#' @family decoding functions
#'
#' @param delta initial distribution of length N, or matrix of dimension c(k,N) for k independent tracks, if \code{trackID} is provided
#' 
#' This could e.g. be the periodically stationary distribution (for each track).
#' @param Gamma array of transition probability matrices for each time point in the cycle of dimension c(N,N,L), where L is the length of the cycle
#' @param allprobs matrix of state-dependent probabilities/ density values of dimension c(n, N)
#' @param tod (Integer valued) variable for cycle indexing in 1, ..., L, mapping the data index to a generalised time of day (length n)
#' 
#' For half-hourly data L = 48. It could, however, also be day of year for daily data and L = 365.
#' @param trackID optional vector of k track IDs, if multiple tracks need to be decoded separately
#' @param mod optional model object containing initial distribution \code{delta}, transition probability matrix \code{Gamma}, matrix of state-dependent probabilities \code{allprobs}, and potentially a \code{trackID} variable
#' 
#' If you are using automatic differentiation either with \code{RTMB::MakeADFun} or \code{\link{qreml}} and include \code{\link{forward_p}} in your likelihood function, the objects needed for state decoding are automatically reported after model fitting.
#' Hence, you can pass the model object obtained from running \code{RTMB::report()} or from \code{\link{qreml}} directly to this function.
#'
#' @return vector of decoded states of length n
#' @export
#'
#' @examples
#' delta = c(0.5, 0.5)
#' beta = matrix(c(-2, 1, -1,
#'                 -2, -1, 1), nrow = 2, byrow = TRUE)
#' Gamma = tpm_p(1:24, 24, beta)
#' 
#' tod = rep(1:24, 5)
#' n = length(tod)
#' 
#' allprobs = matrix(runif(2*n), nrow = n, ncol = 2)
#' states = viterbi_p(delta, Gamma, allprobs, tod)
viterbi_p = function(delta, Gamma, allprobs, tod, trackID = NULL,
                     mod = NULL){
  
  # check if a model with delta, Gamma and allprobs is provided
  if(!is.null(mod)){
    if(is.null(mod$delta)){
      stop("Model object contains no initial distribution.")
    }
    if(is.null(mod$Gamma)){
      stop("Model object contains no transition matrix.")
    }
    if(is.null(mod$allprobs)){
      stop("Model object contains no state-dependent probabilities.")
    }
    if(is.null(mod$tod)){
      stop("Model object contains no cyclic indexing variable.")
    }
    
    # if suitable model object is provided, overwrite inputs with model object
    delta = mod$delta
    Gamma = mod$Gamma
    allprobs = mod$allprobs
    tod = mod$tod
    trackID = mod$trackID
  }
  
  n = nrow(allprobs)
  N = ncol(allprobs)
  
  # creating repeating Gamma array from L unique tpms
  Gammanew = Gamma[,,tod]
  
  if(is.null(trackID)){
    Gammanew = Gammanew[,,-1]
  }
  
  # as.integer(viterbi_g_cpp(allprobs, delta, Gammanew))  
  viterbi_g(delta, Gammanew, allprobs, trackID) 
}



# Local decoding ----------------------------------------------------------


#' Calculate conditional local state probabilities for homogeneous HMMs
#' 
#' Computes
#' \deqn{\Pr(S_t = j \mid X_1, ..., X_T)}
#' for homogeneous HMMs
#'
#' @family decoding functions
#' 
#' @param delta initial or stationary distribution of length N, or matrix of dimension c(k,N) for k independent tracks, if \code{trackID} is provided
#' @param Gamma transition probability matrix of dimension c(N,N), or array of k transition probability matrices of dimension c(N,N,k), if \code{trackID} is provided
#' @param allprobs matrix of state-dependent probabilities/ density values of dimension c(n, N)
#' @param trackID optional vector of length n containing IDs
#' 
#' If provided, the total log-likelihood will be the sum of each track's likelihood contribution.
#' In this case, \code{Gamma} can be a matrix, leading to the same transition probabilities for each track, or an array of dimension c(N,N,k), with one (homogeneous) transition probability matrix for each track.
#' Furthermore, instead of a single vector \code{delta} corresponding to the initial distribution, a \code{delta} matrix of initial distributions, of dimension c(k,N), can be provided, such that each track starts with it's own initial distribution.
#' @param mod optional model object containing initial distribution \code{delta}, transition probability matrix \code{Gamma}, matrix of state-dependent probabilities \code{allprobs}, and potentially a \code{trackID} variable
#' 
#' If you are using automatic differentiation either with \code{RTMB::MakeADFun} or \code{\link{qreml}} and include \code{\link{forward}} in your likelihood function, the objects needed for state decoding are automatically reported after model fitting.
#' Hence, you can pass the model object obtained from running \code{RTMB::report()} or from \code{\link{qreml}} directly to this function.
#'
#' @return matrix of conditional state probabilities of dimension c(n,N)
#' @export
#'
#' @examples
#' Gamma = tpm(c(-1,-2))
#' delta = stationary(Gamma)
#' allprobs = matrix(runif(10), nrow = 10, ncol = 2)
#' 
#' probs = stateprobs(delta, Gamma, allprobs)
stateprobs = function(delta, Gamma, allprobs, trackID = NULL,
                      mod = NULL) {
  
  # check if a model with delta, Gamma and allprobs is provided
  if(!is.null(mod)){
    if(is.null(mod$delta)){
      stop("Model object contains no initial distribution.")
    }
    if(is.null(mod$Gamma)){
      stop("Model object contains no transition matrix.")
    }
    if(is.null(mod$allprobs)){
      stop("Model object contains no state-dependent probabilities.")
    }
    
    # if suitable model object is provided, overwrite inputs with model object
    delta = mod$delta
    Gamma = mod$Gamma
    allprobs = mod$allprobs
    trackID = mod$trackID
  }
  
  n = nrow(allprobs)
  N = ncol(allprobs)
  
  # inflating Gamma to use stateprobs_g
  if(is.null(trackID)){
    Gammanew = array(Gamma, dim = c(N, N, n-1))
  } else{
    uID = unique(trackID)
    k = length(uID) # number of tracks
    
    if(is.matrix(Gamma)){
      Gammanew = array(Gamma, dim = c(N, N, n))
    } else if(is.array(Gamma)){
      if(dim(Gamma)[3] != k) stop("Number of distinct transition matrices does not match the number of tracks.")
      
      ## construct integer trackID
      integerID = match(trackID, uID)
      
      Gammanew = Gamma[,,integerID]
      Gammanew = Gammanew[,, -1] # ignoring first slice
    }
  }
  
  stateprobs_g(delta, Gammanew, allprobs, trackID)
}


#' Calculate conditional local state probabilities for inhomogeneous HMMs
#' 
#' Computes
#' \deqn{\Pr(S_t = j \mid X_1, ..., X_T)}
#' for inhomogeneous HMMs
#' 
#' @family decoding functions
#' 
#' @param delta initial or stationary distribution of length N, or matrix of dimension c(k,N) for k independent tracks, if \code{trackID} is provided
#' @param Gamma array of transition probability matrices of dimension c(N,N,n-1), as in a time series of length n, there are only n-1 transitions
#' 
#' If an array of dimension c(N,N,n) for a single track is provided, the first slice will be ignored.
#' 
#' If \code{trackID} is provided, \code{Gamma} needs to be an array of dimension c(N,N,n), where n is the number of rows in \code{allprobs}. Then for each track the first transition matrix will be ignored.
#' @param allprobs matrix of state-dependent probabilities/ density values of dimension c(n, N)
#' @param trackID optional vector of k track IDs, if multiple tracks need to be decoded separately
#' @param mod optional model object containing initial distribution \code{delta}, transition probability matrix \code{Gamma}, matrix of state-dependent probabilities \code{allprobs}, and potentially a \code{trackID} variable
#' 
#' If you are using automatic differentiation either with \code{RTMB::MakeADFun} or \code{\link{qreml}} and include \code{\link{forward_g}} in your likelihood function, the objects needed for state decoding are automatically reported after model fitting.
#' Hence, you can pass the model object obtained from running \code{RTMB::report()} or from \code{\link{qreml}} directly to this function.
#'
#' @return matrix of conditional state probabilities of dimension c(n,N)
#' @export
#'
#' @examples
#' Gamma = tpm_g(runif(10), matrix(c(-1,-1,1,-2), nrow = 2, byrow = TRUE))
#' delta = c(0.5, 0.5)
#' allprobs = matrix(runif(20), nrow = 10, ncol = 2)
#' 
#' probs = stateprobs_g(delta, Gamma[,,-1], allprobs)
stateprobs_g = function(delta, Gamma, allprobs, trackID = NULL,
                        mod = NULL) {
  
  # check if a model with delta, Gamma and allprobs is provided
  if(!is.null(mod)){
    if(is.null(mod$delta)){
      stop("Model object contains no initial distribution.")
    }
    if(is.null(mod$Gamma)){
      stop("Model object contains no transition matrix.")
    }
    if(is.null(mod$allprobs)){
      stop("Model object contains no state-dependent probabilities.")
    }
    
    # if suitable model object is provided, overwrite inputs with model object
    delta = mod$delta
    # Gamma = mod$Gamma
    Gamma = mod$Gamma[,, -1] # ignoring first slice
    allprobs = mod$allprobs
    trackID = mod$trackID
  }
  
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
    for(i in seq_len(length(uID))){
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
        stateprobs[id_i, ] = exp(lalpha + lbeta - llk)
        
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
    
    c = max(lalpha[nrow(lalpha),])
    llk = c + log(sum(exp(lalpha[nrow(lalpha),] - c)))
    
    stateprobs = exp(lalpha + lbeta - llk)
  }
  
  # rowSums should already be one, but just to be safe
  stateprobs <- stateprobs / rowSums(stateprobs)
  
  colnames(stateprobs) <- paste0("S", 1:N)
  stateprobs
}


#' Calculate conditional local state probabilities for periodically inhomogeneous HMMs
#' 
#' Computes
#' \deqn{\Pr(S_t = j \mid X_1, ..., X_T)}
#' for periodically inhomogeneous HMMs
#'  
#' @family decoding functions
#' 
#' @param delta initial or stationary distribution of length N, or matrix of dimension c(k,N) for k independent tracks, if \code{trackID} is provided
#'
#' This could e.g. be the periodically stationary distribution (for each track) as computed by \code{\link{stationary_p}}.
#' @param Gamma array of transition probability matrices for each time point in the cycle of dimension c(N,N,L), where L is the length of the cycle.
#' @param allprobs matrix of state-dependent probabilities/ density values of dimension c(n, N)
#' @param tod (Integer valued) variable for cycle indexing in 1, ..., L, mapping the data index to a generalised time of day (length n).
#' For half-hourly data L = 48. It could, however, also be day of year for daily data and L = 365.
#' @param trackID optional vector of k track IDs, if multiple tracks need to be decoded separately
#' @param mod optional model object containing initial distribution \code{delta}, transition probability matrix \code{Gamma}, matrix of state-dependent probabilities \code{allprobs}, and potentially a \code{trackID} variable
#' 
#' If you are using automatic differentiation either with \code{RTMB::MakeADFun} or \code{\link{qreml}} and include \code{\link{forward_p}} in your likelihood function, the objects needed for state decoding are automatically reported after model fitting.
#' Hence, you can pass the model object obtained from running \code{RTMB::report()} or from \code{\link{qreml}} directly to this function.
#' 
#' @return matrix of conditional state probabilities of dimension c(n,N)
#' @export
#'
#' @examples
#' L = 24
#' beta = matrix(c(-1, 2, -1, -2, 1, -1), nrow = 2, byrow = TRUE)
#' Gamma = tpm_p(1:L, L, beta, degree = 1)
#' delta = stationary_p(Gamma, 1)
#' allprobs = matrix(runif(200), nrow = 100, ncol = 2)
#' tod = rep(1:24, 5)[1:100]
#' 
#' probs = stateprobs_p(delta, Gamma, allprobs, tod)

stateprobs_p = function(delta, Gamma, allprobs, tod, trackID = NULL,
                        mod = NULL) {
  
  # check if a model with delta, Gamma and allprobs is provided
  if(!is.null(mod)){
    if(is.null(mod$delta)){
      stop("Model object contains no initial distribution.")
    }
    if(is.null(mod$Gamma)){
      stop("Model object contains no transition matrix.")
    }
    if(is.null(mod$allprobs)){
      stop("Model object contains no state-dependent probabilities.")
    }
    if(is.null(mod$tod)){
      stop("Model object contains no cyclic indexing variable.")
    }
    
    # if suitable model object is provided, overwrite inputs with model object
    delta = mod$delta
    Gamma = mod$Gamma
    allprobs = mod$allprobs
    tod = mod$tod
    trackID = mod$trackID
  }
  
  n = nrow(allprobs)
  N = ncol(allprobs)
  
  Gammanew = Gamma[,,tod] # select the transition matrix for the current time of day
  
  if(is.null(trackID)){
    Gammanew = Gammanew[,,-1]
  }
  
  stateprobs_g(delta, Gammanew, allprobs, trackID)
}

