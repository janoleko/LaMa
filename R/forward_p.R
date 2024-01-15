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
#' ## generating data from periodic 2-state HMM
#' mu = c(0, 6)
#' sigma = c(2, 4)
#' beta = matrix(c(-2,-2,1,-1, 1, -1),nrow=2)
#' delta = c(0.5, 0.5)
#' # simulation
#' n = 2000
#' s = x = rep(NA, n)
#' tod = rep(1:24, ceiling(2000/24))
#' s[1] = sample(1:2, 1, prob = delta)
#' x[1] = stats::rnorm(1, mu[s[1]], sigma[s[1]])
#' # 24 unique t.p.m.s
#' Gamma = array(dim = c(2,2,24))
#' for(t in 1:24){
#'   G = diag(2)
#'   G[!G] = exp(beta[,1]+beta[,2]*sin(2*pi*t/24)+
#'     beta[,3]*cos(2*pi*t/24)) # trigonometric link
#'   Gamma[,,t] = G / rowSums(G)
#' }
#' for(t in 2:n){
#'   s[t] = sample(1:2, 1, prob = Gamma[s[t-1],,tod[t]])
#'   x[t] = stats::rnorm(1, mu[s[t]], sigma[s[t]])
#' }
#' 
#' ## negative log likelihood function
#' mllk = function(theta.star, x, tod){
#'   # parameter transformations for unconstraint optimization
#'   beta = matrix(theta.star[1:6], 2, 3)
#'   Gamma = array(dim=c(2,2,24))
#'   for(t in 1:24){
#'     G = diag(2)
#'     G[!G] = exp(beta[,1]+beta[,2]*sin(2*pi*t/24)+
#'     beta[,3]*cos(2*pi*t/24)) # trigonometric link
#'     Gamma[,,t] = G / rowSums(G)
#'   }
#'   delta = c(plogis(theta.star[7]), 1-plogis(theta.star[5]))
#'   mu = theta.star[8:9]
#'   sigma = exp(theta.star[10:11])
#'   # calculate all state-dependent probabilities
#'   allprobs = matrix(1, length(x), 2)
#'   for(j in 1:2){ allprobs[,j] = stats::dnorm(x, mu[j], sigma[j]) }
#'   # return negative for minimization
#'   -Lcpp::forward_p(delta, Gamma, allprobs, tod)
#' }
#' 
#' ## fitting an HMM to the data
#' theta.star = c(-2,-2,1,-1,1,-1,0,0,5,log(2),log(3))
#' mod = stats::nlm(mllk, theta.star, x = x, tod = tod)
#'
forward_p = function(delta, Gamma, allprobs, tod){
  if(min(tod)==1){
    tod = tod-1
  } 
  forward_cpp_p(allprobs, delta, Gamma, tod)
}