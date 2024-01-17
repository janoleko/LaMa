#' Trigonometric Basis Expansion
#'
#' Given a periodically varying variable such as time of day or day of year and the associated cycle length, this function performs a basis expansion to efficiently calculate a linear predictor of the form \cr \cr
#' \eqn{ 
#'  \eta = \beta_0 + \sum_{k=1}^K \bigl( \beta_{1k} \sin(\frac{2 \pi k t}{L}) + \beta_{2k} \cos(\frac{2 \pi k t}{L}) \bigr). 
#'  } \cr \cr
#'  This is relevant for modeling e.g. diurnal variation and the flexibility can be increased by adding smaller frequencies (i.e. increasing \eqn{K}).
#'  
#' @param tod Time variable, describing the time point in a cycle. Could for example be time of day (between 0 and 24) or day of year.
#' @param L Length of one cycle on the scale of the time variable. For time of day, this would be 24.
#' @param degree Degree K of the trigonometric link above. Increasing K increases the flexibility.
#'
#' @return A design matrix (without intercept column of ones), ordered as sin1, cos1, sin2, cos2, ...
#' @export
#'
#' @examples
#' ## hourly data
#' tod = rep(1:24, 10)
#' Z = trigBasisExp(tod, L = 24, degree = 2)
#' 
#' ## half-hourly data
#' tod = rep(1:48/2, 10) # in [0,24] -> L = 24
#' Z1 = trigBasisExp(tod, L = 24, degree = 3)
#' 
#' tod = rep(1:48, 10) # in [1,48] -> L = 48
#' Z2 = trigBasisExp(tod, L = 48, degree = 3)
#' 
#' Z1 - Z2
#' # The latter two are equivalent specifications!

trigBasisExp = function(tod, L = 24, degree = 1){
  n = length(tod)
  Z = matrix(nrow = n, ncol = 2*degree)
  inner = 2*pi*tod/L
  for(k in 1:degree){
    Z[,2*(k-1)+1:2] = cbind(sin(inner*k), cos(inner*k))
  }
  colnames(Z) = paste0(c("sin_", "cos_"), rep(1:degree, each = 2))
  Z
}
