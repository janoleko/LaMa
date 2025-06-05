#' AD-compatible minimum and maximum functions
#' 
#' These functions compute the parallel minimum/ maximum of two vector-valued inputs and are compatible with automatic differentiation using \code{RTMB}.
#'
#' @param x first vector
#' @param y second vector
#'
#' @returns \code{min2} returns the parallel minimum and \code{max2} the parallel maximum of \code{x} and \code{y}
#'
#' @examples
#' x <- c(1, 4, 8, 2)
#' y <- c(2, 5, 3, 7)
#' min2(x, y)
#' max2(x, y)
#' @name minmax
NULL

#' @rdname minmax
#' @export
min2 <- function(x,y){
  (x + y - abs(x - y)) / 2
}
#' @rdname minmax
#' @export
max2 = function(x,y){
  (x + y + abs(x - y)) / 2
}

#' Calculate the index of the first observation of each track based on an ID variable
#' 
#' Function to conveniently calculate the trackInd variable that is needed internally when fitting a model to longitudinal data with multiple tracks.
#'
#' @param ID ID variable of track IDs that is of the same length as the data to be analysed
#'
#' @return A vector of indices of the first observation of each track which can be passed to the forward and forward_g to sum likelihood contributions of each track
#' @export
#'
#' @examples
#' uniqueID = c("Animal1", "Animal2", "Animal3")
#' ID = rep(uniqueID, c(100, 200, 300))
#' trackInd = calc_trackInd(ID)
calc_trackInd = function(ID){
  
  if(is.factor(ID)){
    ID = as.vector(ID)
  }
  
  if(!is.vector(ID)){
    stop("ID must be a vector")
  }
  
  RLE = rle(ID)
  trackInd = c(1, cumsum(RLE$lengths[-length(RLE$lengths)]) + 1)
  
  return(trackInd)
}

# helper functions for penalty and qreml
reshape_lambda <- function(num_elements, lambda) {
  start <- 1
  result <- lapply(num_elements, function(len) {
    # Extract sub-vector from lambda based on the number of elements
    sub_vector <- lambda[start:(start + len - 1)]
    start <<- start + len
    return(sub_vector)
  })
  return(result)
}
map_lambda = function(lambda, map){
  lambda_mapped = numeric(length(levels(map)))
  lambda_levels = levels(map)
  for(l in seq_along(lambda_levels)){
    thislambda <- lambda[which(map == lambda_levels[l])]
    if(is.numeric(thislambda)){ # if numeric, take mean
      lambda_mapped[l] = mean(thislambda)
    } else if(is.character(thislambda)){ # if character, collapse with &
      lambda_mapped[l] = paste(thislambda, collapse = "&")
    }
  }
  names(lambda_mapped) = levels(map)
  return(lambda_mapped)
}
unmap_lambda = function(lambda_mapped, map, lambda0){
  lambda = rep(NA, length(lambda0))
  lambda_levels = levels(map)
  for(l in seq_along(lambda_levels)){
    lambda[which(map == lambda_levels[l])] = lambda_mapped[l]
  }
  NA_ind = which(is.na(lambda))
  lambda[NA_ind] = lambda0[NA_ind]
  return(lambda)
}

#' Computes generalised determinant
#'
#' @param x square matrix
#' @param eps eigenvalues smaller than this will be treated as zero
#' @param log logical. If \code{TRUE}, the log-determinant is returned. If \code{FALSE}, the determinant is returned.
#'
#' @return generalised log-determinant of \code{x}
#' 
#' @importFrom RTMB eigen
gdeterminant <- function(x, eps = NULL, log = TRUE) {
  if(is.null(x)) {
    return(NULL)
  } else {
    # Compute sum of log of non-zero eigenvalues
    # (i.e., log generalized determinant)
    eigenpairs <- eigen(x)
    eigenvalues <- eigenpairs$values
    # if no tolerance provided: largest eigenvalue times machine precision
    if(is.null(eps)) {
      eps <- max(eigenvalues) * .Machine$double.eps
    }
    logdet <- sum(log(eigenvalues[eigenvalues > eps]))
    if(!log) {
      return(exp(logdet))
    } else{
      return(logdet)
    }        
  }
}


logspace_add <- function(x){
  # Computes the log of the sum of exponentials of the input vector x
  # This is useful for numerical stability when dealing with large numbers
  max_x <- max(x)
  max_x + log(sum(exp(x - max_x)))
}
