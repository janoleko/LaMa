#' Forward algorithm with (only) periodically varying transition probability matrix
#'
#' When the transition probability matrix only varies periodically (e.g. as a function of time of day), there are only L unique matrices if L is the period length (e.g. L = 24 for hourly data and time-of-day variation).
#' Thus it is much more efficient to only calculate these L matrices and index them by time of day instead of calculating such a matrix for each index in the data set.
#' This function allows for exactly this, by only expecting a Gamma matrix for each time point in a day, and an integer valued (1, ..., L) time of day variable that maps the data index to the according time of day.
#'
#' @param delta Initial or periodically stationary distribution (of length N)
#' @param Gamma Pre-calculated array of periodic Gamma matrices (of dimension c(N,N,L))
#' @param allprobs allprobs matrix (of dimension c(n, N))
#' @param tod (Integer valued) time variable in 1, ..., L, mapping the data index to a generalized time of day (length n).
#' For half-hourly data L = 48. It could however also be day of year when L = 365.
#'
#' @return Log-likelihood for given data and parameters
#' @export
#'
#' @examples
forward_p = function(delta, Gamma, allprobs, tod){
  forward_cpp_p(allprobs, delta, Gamma, tod)
}