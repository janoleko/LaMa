# differentiable max function
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

# helper function for penalty and qreml
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

# helper function to compute generalized determinant
gdeterminant <- function(x, 
                         eps = 1e-10, # eigenvalues smaller than this will be treated as zero
                         log = TRUE){
  svd = eigen(x)
  values = svd$values
  
  logdet = sum(log(values[values > eps]))
  
  if(!log){
    return(exp(logdet))
  } else{
    return(logdet)
  }
}