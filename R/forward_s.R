#' \href{https://www.taylorfrancis.com/books/mono/10.1201/b20790/hidden-markov-models-time-series-walter-zucchini-iain-macdonald-roland-langrock}{Forward algorithm} for hidden semi-Markov models with homogeneous transition probability matrix
#'
#' Hidden semi-Markov models (HSMMs) are a flexible extension of HMMs. 
#' For direct numerical maximum likelhood estimation, HSMMs can be represented as HMMs on an enlarged state space (of size \eqn{M}) and with structured transition probabilities.
#'
#' @param delta Initial or stationary distribution of length M (where M is the number of approximating states)
#' @param Gamma Transition probability matrix of dimension c(M,M)
#' @param allprobs Matrix of state-dependent probabilities/ density values of dimension c(n, N), where N is the number of semi-Markovian states.
#' This will automatically be converted to the appropriate dimension.
#' @param sizes State aggregate sizes that are used for the approximation of the semi-Markov chain.
#'
#' @return Log-likelihood for given data and parameters
#' @export
#'
#' @examples
#' ## generating data from homogeneous 2-state HSMM
#' mu = c(0, 6)
#' lambda = c(6, 12)
#' omega = matrix(c(0,1,1,0), nrow = 2, byrow = TRUE)
#' # simulation
#' # for a 2-state HSMM the embedded chain always alternates between 1 and 2
#' s = rep(1:2, 100)
#' C = x = numeric(0)
#' for(t in 1:100){
#'   dt = rpois(1, lambda[s[t]])+1 # shifted Poisson
#'   C = c(C, rep(s[t], dt))
#'   x = c(x, rnorm(dt, mu[s[t]], 1.5)) # fixed sd 2 for both states
#' }
#' 
#' ## negative log likelihood function
#' mllk = function(theta.star, x, sizes){
#'   # parameter transformations for unconstraint optimization
#'   omega = matrix(c(0,1,1,0), nrow = 2, byrow = TRUE) # omega fixed (2-states)
#'   lambda = exp(theta.star[1:2]) # dwell time means
#'   dm = list(dpois(1:sizes[1]-1, lambda[1]), dpois(1:sizes[2]-1, lambda[2]))
#'   Gamma = Lcpp::tpm_hsmm(omega, dm)
#'   delta = Lcpp::stationary(Gamma) # stationary
#'   mu = theta.star[3:4]
#'   sigma = exp(theta.star[5:6])
#'   # calculate all state-dependent probabilities
#'   allprobs = matrix(1, length(x), 2)
#'   for(j in 1:2){ allprobs[,j] = stats::dnorm(x, mu[j], sigma[j]) }
#'   # return negative for minimization
#'   -Lcpp::forward_s(delta, Gamma, allprobs, sizes)
#' }
#' 
#' ## fitting an HSMM to the data
#' theta.star = c(log(5), log(10), 1, 4, log(2), log(2))
#' mod = stats::nlm(mllk, theta.star, x = x, sizes = c(20, 30), stepmax = 5)
forward_s = function(delta, Gamma, allprobs, sizes){
  forward_cpp_s(allprobs, delta, Gamma, sizes)
}