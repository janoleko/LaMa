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
#' @param x symmetric matrix
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

## log(sum(exp(x))) # with numerical stability
logspace_add <- function(x){
  # Computes the log of the sum of exponentials of the input vector x
  # This is useful for numerical stability when dealing with large numbers
  max_x <- max(x)
  max_x + log(sum(exp(x - max_x)))
}


## function to construct unimodality constraint matrix (list)
construct_C <- function(m, # beta mode index (vector of length N)
                        k, # basis dimension
                        exclude_last = TRUE)
{
  # construct constraint matrices
  m <- as.numeric(m)
  N <- length(m) # states
  if(any(m > k)) stop("'m' must be less than or equal to 'k'")
  D <- diff(diag(k)) # first-order difference matrix
  C <- list() # initialise list
  for(i in 1:N) {
    Ci <- D # initialise with first-order difference matrix
    if(m[i] < k){ # if there is a block that needs to be flipped, do that
      Ci[m[i]:(k-1), ] <- -Ci[m[i]:(k-1),]
    }
    if(exclude_last) Ci <- Ci[,-k] # maybe exclude last column because multiplied by zero
    C[[i]] <- Ci
  }
  return(C)
}

#' Smooth approximations to max(x, 0) and min(x, 0)
#'
#' @param x a vector of values
#' @param rho smoothing parameter, larger values lead to closer approximation
#'
#' @return the approximate maximum or minimum of x and 0
#'
#' @examples
#' x <- seq(-1, 1, by = 0.1)
#' min0_smooth(x)
#' max0_smooth(x)
#' @name minmax0_smooth
NULL

#' @rdname minmax0_smooth
#' @export
max0_smooth <- function(x, rho = 20){
  (1 / rho) * log1p(exp(rho * x))
}

#' @rdname minmax0_smooth
#' @export
min0_smooth <- function(x, rho = 20){
  # (-1 / rho) * log(1 + exp(-rho * x))
  (-1 / rho) * log1p(exp(-rho * x))
}


#' Sparsity-retaining matrix multiplication
#' 
#' Currently, standard matrix multiplication destroys automatic sparsity detection by \code{RTMB} which is essential for models with high-dimensional random effects.
#' This version retains said sparsity. It may be slightly slower when constructing the AD tape.
#'
#' @param A matrix of dimension n x p
#' @param B matrix of dimension p x m
#'
#' @return the matrix product of A and B, which is of dimension n x m
#' @export
#'
#' @examples
#' A <- matrix(1:6, nrow = 2, ncol = 3)
#' B <- matrix(7:12, nrow = 3, ncol = 2)
#' A %sp% B
`%sp%` <- function(A, B) {
  A <- as.matrix(A)
  B <- as.matrix(B)
  
  # n <- nrow(A)
  p <- ncol(A)
  q <- nrow(B)
  m <- ncol(B)
  
  if (p != q) {
    stop("Non-conformable matrices: ncol(A) must equal nrow(B)")
  }
  
  At <- t(A)  # p x n
  sapply(1:m, function(j) colSums(At * B[, j])) # sums of element-wise multiplication retain sparsity
}
