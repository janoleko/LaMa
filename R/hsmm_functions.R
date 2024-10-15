#' \href{https://www.taylorfrancis.com/books/mono/10.1201/b20790/hidden-markov-models-time-series-walter-zucchini-iain-macdonald-roland-langrock}{Forward algorithm} for homogeneous hidden semi-Markov models
#'
#' @description
#' Hidden semi-Markov models (HSMMs) are a flexible extension of HMMs, where the state duration distribution is explicitly modelled by a distribution on the positive integers.
#' For direct numerical maximum likelhood estimation, HSMMs can be represented as HMMs on an enlarged state space (of size \eqn{M}) and with structured transition probabilities.
#' 
#' This function is designed to be used with automatic differentiation based on the \code{R} package \code{RTMB}. It will be very slow without it.
#'
#' @param dm List of length N containing vectors of dwell-time probability mass functions (PMFs) for each state. The vector lengths correspond to the approximating state aggregate sizes, hence there should be little probablity mass not covered by these.
#' @param omega Matrix of dimension c(N,N) of conditional transition probabilites. This contains the transition probabilities given the current state is left. Hence, the diagonal elements need to be zero and the rows need to sum to one.
#' @param allprobs Matrix of state-dependent probabilities/ density values of dimension c(n, N)
#' @param trackID Optional vector of length n containing k unique IDs. If provided, the total log-likelihood will be the sum of each track's likelihood contribution.
#' In this case, \code{dm} can be a nested list, where the top layer contains k \code{dm} lists as described above. \code{omega} can then also be an array of dimension c(N,N,k) with one conditional transition probability matrix for each track.
#' Furthermore, instead of a single vector \code{delta} corresponding to the initial distribution, a delta matrix of initial distributions, of dimension c(k,N), can be provided, such that each track starts with it's own initial distribution.
#' @param delta Optional vector of initial state probabilities of length N. By default, the stationary distribution is computed (which is typically recommended).
#' @param eps Small value to avoid numerical issues in the approximating transition matrix construction.
#' @param report Logical, indicating whether initial distribution, approximating t.p.m. and allprobs matrix should be reported from the fitted model. Defaults to TRUE.
#'
#' @return Log-likelihood for given data and parameters
#' @export
#' @import RTMB
#'
#' @examples
#' # currently no examples
forward_hsmm <- function(dm, omega, allprobs,
                         trackID = NULL, delta = NULL, eps = 1e-20, report = TRUE){
  # overloading assignment operators, currently necessary
  "[<-" <- ADoverload("[<-")
  "c" <- ADoverload("c")
  "diag<-" <- ADoverload("diag<-")
  ################################
  
  agsizes = sapply(dm, length)
  
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
#' @description
#' Hidden semi-Markov models (HSMMs) are a flexible extension of HMMs, where the state duration distribution is explicitly modelled by a distribution on the positive integers. This function can be used to fit HSMMs where the state-duration distribution and/ or the conditional transition probabilities vary with covariates.
#' For direct numerical maximum likelhood estimation, HSMMs can be represented as HMMs on an enlarged state space (of size \eqn{M}) and with structured transition probabilities.
#' 
#' This function is designed to be used with automatic differentiation based on the \code{R} package \code{RTMB}. It will be very slow without it.
#'
#' @param dm List of length N containing matrices (or vectors) of dwell-time probability mass functions (PMFs) for each state.
#' 
#' If the dwell-time PMFs are constant, the vectors are the PMF of the dwell-time distribution fixed in time. The vector lengths correspond to the approximating state aggregate sizes, hence there should be little probablity mass not covered by these.
#' 
#' If the dwell-time PMFs are inhomogeneous, the matrices need to have n rows, where n is the number of observations. The number of columns again correponds to the size of the approximating state aggregates.
#' 
#' In the latter case, the first \code{max(sapply(dm, ncol)) - 1} observations will not be used because the first t.p.m. needs to be computed based on the first \code{max(sapply(dm, ncol))} covariate values (represented by \code{dm}).
#' @param omega matrix of dimension c(N,N) or array of dimension c(N,N,n) of conditional transition probabilites. This contains the transition probabilities given the current state is left. Hence, the diagonal elements need to be zero and the rows need to sum to one.
#' @param allprobs Matrix of state-dependent probabilities/ density values of dimension c(n, N)
#' @param trackID Optional vector of length n containing k unique IDs. If provided, the total log-likelihood will be the sum of each track's likelihood contribution.
#' Instead of a single vector \code{delta} corresponding to the initial distribution, a \code{delta} matrix of initial distributions, of dimension c(k,N), can be provided, such that each track starts with it's own initial distribution.
#' @param delta Optional vector of initial state probabilities of length N. By default, instead of this, the stationary distribution is computed corresponding to the first approximating t.p.m. of each track is computed. Contrary to the homogeneous case, this is not theoretically motivated but just for convenience.
#' @param eps Small value to avoid numerical issues in the approximating transition matrix construction.
#' @param report Logical, indicating whether initial distribution, approximating t.p.m. and allprobs matrix should be reported from the fitted model. Defaults to TRUE.
#'
#' @return Log-likelihood for given data and parameters
#' @export
#' @import RTMB
#'
#' @examples
#' # currently no examples
forward_ihsmm <- function(dm, omega, allprobs,
                          trackID = NULL, delta = NULL, eps = 1e-20, report = TRUE){
  
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
  
  stationary = is.null(delta) # if delta is not provided, stationary distribution needs to be computed
  
  if(is.null(trackID)) {
    
    ## first case: dm varies -> truely inhomogeneous HSMM
    if(is.matrix(dm[[1]])){
      if(nrow(dm[[1]]) != n) stop("dm needs to be a list of length N, either containing vectors or matrices of dwell-time PMFs.")
      
      ## compute approximating tpm
      Gamma_sparse = tpm_ihsmm(omega, dm, eps = eps)
      
      ## if stationary, compute initial stationary distribution
      if(stationary){
        delta = stationary_sparse(Gamma_sparse[[1]])
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
      foo = delta_sparse %*% diag(rep(allprobs[maxag,], times = agsizes))
      sumfoo = sum(foo)
      phi = foo / sumfoo
      l = log(sumfoo)
      
      for(t in (maxag + 1):nrow(allprobs)) {
        # foo = phi %*% Gamma_sparse %*% Matrix::Diagonal(x = rep(allprobs[t,], times = agsizes))
        foo = phi %*% Gamma_sparse[[t-maxag]] %*% diag(rep(allprobs[t,], times = agsizes))
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
        foo = delta_i %*% diag(rep(allprobs[ind[maxag],], times = agsizes))
        sumfoo = sum(foo)
        phi = foo / sumfoo
        l_i = log(sumfoo)
        
        for(t in (maxag + 1):length(ind)) {
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
#' @description
#' Hidden semi-Markov models (HSMMs) are a flexible extension of HMMs, where the state duration distribution is explicitly modelled by a distribution on the positive integers.
#'
#' This function can be used to fit HSMMs where the state-duration distribution and/ or the conditional transition probabilities vary periodically.
#' In the special case of periodic variation (as compared to arbitrary covariate influence), this version is to be preferred over \code{\link{forward_ihsmm}} because it computes the correct periodically stationary distribution and no observations are lost for the approximation.
#' For direct numerical maximum likelhood estimation, HSMMs can be represented as HMMs on an enlarged state space (of size \eqn{M}) and with structured transition probabilities.
#' 
#' This function is designed to be used with automatic differentiation based on the \code{R} package \code{RTMB}. It will be very slow without it.
#'
#' @param dm List of length N containing matrices (or vectors) of dwell-time probability mass functions (PMFs) for each state.
#'
#' If the dwell-time PMFs are constant, the vectors are the PMF of the dwell-time distribution fixed in time. The vector lengths correspond to the approximating state aggregate sizes, hence there should be little probablity mass not covered by these.
#' 
#' If the dwell-time PMFs are inhomogeneous, the matrices need to have L rows, where L is the cycle length. The number of columns again correpond to the size of the approximating state aggregates.
#' @param omega Matrix of dimension c(N,N) or array of dimension c(N,N,L) of conditional transition probabilites. This contains the transition probabilities given the current state is left. Hence, the diagonal elements need to be zero and the rows need to sum to one.
#' @param allprobs Matrix of state-dependent probabilities/ density values of dimension c(n, N)
#' @param tod (Integer valued) time variable in 1, ..., L, mapping the data index to a generalized time of day (length n).
#' For half-hourly data L = 48. It could, however, also be day of year for daily data and L = 365.
#' @param trackID Optional vector of length n containing k unique IDs. If provided, the total log-likelihood will be the sum of each track's likelihood contribution.
#' Instead of a single vector \code{delta} corresponding to the initial distribution, a \code{delta} matrix of initial distributions, of dimension c(k,N), can be provided, such that each track starts with it's own initial distribution.
#' @param delta Optional vector of initial state probabilities of length N. By default, instead of this, the stationary distribution is computed corresponding to the first approximating t.p.m. of each track is computed. Contrary to the homogeneous case, this is not theoretically motivated but just for convenience.
#' @param eps Small value to avoid numerical issues in the approximating transition matrix construction.
#' @param report Logical, indicating whether initial distribution, approximating t.p.m. and allprobs matrix should be reported from the fitted model. Defaults to TRUE.
#'
#' @return Log-likelihood for given data and parameters
#' @export
#' @import RTMB
#'
#' @examples
#' # currently no examples
forward_phsmm <- function(dm, omega, allprobs, tod,
                          trackID = NULL, delta = NULL, eps = 1e-20, report = TRUE){
  
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


# differentiable max function
max2 = function(x,y){
  (x + y + abs(x - y)) / 2
}

#' Build the transition probability matrix of an HSMM-approximating HMM
#'
#' @description
#' Hidden semi-Markov models (HSMMs) are a flexible extension of HMMs, where the state duration distribution is explicitly modelled.
#' For direct numerical maximum likelhood estimation, HSMMs can be represented as HMMs on an enlarged state space (of size \eqn{M}) and with structured transition probabilities.
#' This function computes the transition matrix to approximate a given HSMM by an HMM with a larger state space.
#'
#' @param omega Embedded transition probability matrix of dimension c(N,N)
#' @param dm State dwell-time distributions arranged in a list of length(N). Each list element needs to be a vector of length N_i, where N_i is the state aggregate size.
#' @param Fm Optional list of length N containing the cumulative distribution functions of the dwell-time distributions.
#' @param sparse Logical, indicating whether the output should be a sparse matrix. Defaults to TRUE.
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
tpm_hsmm <- function(omega, dm, 
                     Fm = NULL, sparse = TRUE, eps = 1e-20) {
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
    ci = dm[[i]] / (1 - Fm[[i]] + eps)
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
          Gi[1, col_start:(col_start + Nj - 1)] = c(omega[i, j] * ci, rep(0, Nj - 1))
        } else {
          Gi[, col_start] = omega[i, j] * ci
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

#' Build all transition probability matrices of an inhomogeneous-HSMM-approximating HMM
#'
#' Hidden semi-Markov models (HSMMs) are a flexible extension of HMMs. For direct numerical maximum likelhood estimation, HSMMs can be represented as HMMs on an enlarged state space (of size \eqn{M}) and with structured transition probabilities.
#' This function computes the transition matrices of a periodically inhomogeneos HSMMs.
#'
#' @param omega Embedded transition probability matrix
#' Either a matrix of dimension c(N,N) for homogeneous conditional transition probabilities, or an array of dimension c(N,N,n) for inhomogeneous conditional transition probabilities.
#' @param dm State dwell-time distributions arranged in a list of length(N).
#' Each list element needs to be a matrix of dimension c(n, N_i), where each row t is the (approximate) probability mass function of state i at time t.
#' @param eps Rounding value: If an entry of the transition probabily matrix is smaller, than it is rounded to zero.
#'
#' @return A list of dimension length \code{n - max(sapply(dm, ncol))}, containing sparse extended-state-space transition probability matrices for each time point (except the first \code{max(sapply(dm, ncol)) - 1}).
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
                     eps = 1e-20){
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

#' Build all transition probability matrices of an periodic-HSMM-approximating HMM
#'
#' Hidden semi-Markov models (HSMMs) are a flexible extension of HMMs. For direct numerical maximum likelhood estimation, HSMMs can be represented as HMMs on an enlarged state space (of size \eqn{M}) and with structured transition probabilities.
#' This function computes the transition matrices of a periodically inhomogeneos HSMMs.
#'
#' @param omega Embedded transition probability matrix
#' Either a matrix of dimension c(N,N) for homogeneous conditional transition probabilities, or an array of dimension c(N,N,L) for inhomogeneous conditional transition probabilities.
#' @param dm State dwell-time distributions arranged in a list of length(N).
#' Each list element needs to be a matrix of dimension c(L, N_i), where each row t is the (approximate) probability mass function of state i at time t.
#' @param eps Rounding value: If an entry of the transition probabily matrix is smaller, than it is rounded to zero.
#'
#' @return A list of dimension length L, containing sparse extended-state-space transition probability matrices of the approximating HMM for each time point of the cycle.
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

#' Sparse version of \code{\link{stationary}}
#'
#' @description
#' This is function computes the stationary distribution of a Markov chain with a given \strong{sparse} transition probability matrix.
#' Compatible with automatic differentiation by \code{RTMB}
#
#' @param Gamma Sparse transition probability matrix of dimension c(N,N)
#'
#' @return Stationary distribution of the Markov chain with the given transition probability matrix
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
#' @param Gamma List of length L containing sparse transition probability matrices for one cycle.
#' @param t Integer index of the time point in the cycle, for which to calculate the stationary distribution
#' If t is not provided, the function calculates all stationary distributions for each time point in the cycle.
#'
#' @return Either the periodically stationary distribution at time t or all periodically stationary distributions.
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
  "[<-" <- RTMB::ADoverload("[<-")
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

#' Build the embedded transition probability matrix of an HSMM from unconstrained parameter vector
#'
#' @description
#' Hidden semi-Markov models are defined in terms of state durations and an embedded transition probability matrix that contains the conditional transition probabilities given that the current state is left.
#' This function builds the embedded/ conditional transition probability matrix from an unconstraint parameter vector. 
#' For each row of the matrix, the inverse multinomial logistic link is applied.
#' 
#' For a matrix of dimension c(N,N), the number of free off-diagonal elements is N*(N-2), hence also the length of \code{param}.
#'
#' Compatible with automatic differentiation by \code{RTMB}
#'
#' @param param Unconstrained parameter vector of length N*(N-2) where N is the number of states of the Markov chain
#' If the function is called without \code{param}, it will return the conditional transition probability matrix for a 2-state HSMM, which is fixed with 0 diagonal entries and off-diagonal entries equal to 1.
#' @param byrow Logical that indicates if the transition probability matrix should be filled by row. 
#' Defaults to FALSE, but should be set to TRUE if one wants to work with a matrix of beta parameters returned by popular HMM packages like \code{moveHMM}, \code{momentuHMM}, or \code{hmmTMB}.
#'
#' @return Embedded/ conditional transition probability matrix of dimension c(N,N)
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
tpm_emb = function(param = NULL, byrow = FALSE){
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
    
    if(byrow){
      omega = t(omega)
    }
    
    omega = omega / rowSums(omega)
  }
  
  omega
}

#' Build all embedded transition probability matrices of an inhomogeneous HSMM
#'
#' @description
#' Hidden semi-Markov models are defined in terms of state durations and an embedded transition probability matrix that contains the conditional transition probabilities given that the current state is left.
#' This function builds all embedded/ conditional transition probability matrices based on a design and parameter matrix.
#' For each row of the matrix, the inverse multinomial logistic link is applied.
#' 
#' For a matrix of dimension c(N,N), the number of free off-diagonal elements is N*(N-2) which determines the number of rows of the parameter matrix.
#'
#' Compatible with automatic differentiation by \code{RTMB}
#' @param Z Covariate design matrix with or without intercept column, i.e. of dimension c(n, p) or c(n, p+1).
#' If Z has only p columns, an intercept column of ones will be added automatically.
#' @param beta Matrix of coefficients for the off-diagonal elements of the embedded transition probability matrix.
#' Needs to be of dimension c(N*(N-2), p+1), where the first column contains the intercepts.
#' @param byrow Logical that indicates if each transition probability matrix should be filled by row. 
#' Defaults to FALSE, but should be set to TRUE if one wants to work with a matrix of beta parameters returned by popular HMM packages like \code{moveHMM}, \code{momentuHMM}, or \code{hmmTMB}.
#' @param report Logical, indicating whether the coefficient matrix beta should be reported from the fitted model. Defaults to TRUE.
#'
#'
#' @return Array of embedded/ conditional transition probability matrices of dimension c(N,N,n)
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
tpm_emb_g = function(Z, beta, byrow = FALSE, report = TRUE){
  "[<-" <- ADoverload("[<-") # overloading assignment operators, currently necessary
  "c" <- ADoverload("c")
  "diag<-" <- ADoverload("diag<-")
  
  p = ncol(beta) - 1 # number of parameters per state (excluding intercept)
  # for N > 2: N*(N-2) is bijective with solution
  N = as.integer(1 + sqrt(1 + p), 0)
  
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
    
    if(byrow){
      O = t(O)
    }
    
    omega[,,t] = O / rowSums(O)
  }
  
  omega
}
