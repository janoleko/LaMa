#' Calculate the index of the first observation of each track based on an ID variable
#' 
#' Function to conveniently calculate the trackInd variable that is needed when fitting a model to longitudinal data with multiple tracks.
#' 
#' Preferably, this function should not be used inside the likelihood function, as it may slow down the computation speed.
#' Instead, it can be called once and the result can then be passed as an argument to the likelihood function.
#'
#' @param ID ID variable of track IDs that is of the same length as the data to be analyzed
#'
#' @return A vector of indices of the first observation of each track which can be passed to the forward and forward_g to sum likelihood contributions of each track
#' @export
#'
#' @examples
#' uniqueID = c("Animal1", "Animal2", "Animal3")
#' ID = rep(uniqueID, c(100, 200, 300))
#' trackInd = calc_trackInd(ID)
calc_trackInd = function(ID){
  if(!is.vector(ID)){
    stop("ID must be a vector")
  }
  
  RLE = rle(ID)
  trackInd = c(1, cumsum(RLE$lengths[-length(RLE$lengths)]) + 1)

  return(trackInd)
}
