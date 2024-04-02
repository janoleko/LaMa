#' \href{https://www.taylorfrancis.com/books/mono/10.1201/b20790/hidden-markov-models-time-series-walter-zucchini-iain-macdonald-roland-langrock}{Forward algorithm} with homogeneous transition probability matrix
#'
#' @param delta Initial or stationary distribution of length N
#' @param Gamma Transition probability matrix of dimension c(N,N)
#' @param allprobs Matrix of state-dependent probabilities/ density values of dimension c(n, N)
#'
#' @return Log-likelihood for given data and parameters
#' @export
#'
#' @examples
#' ## generating data from homogeneous 2-state HMM
#' mu = c(0, 6)
#' sigma = c(2, 4)
#' Gamma = matrix(c(0.5, 0.05, 0.15, 0.85), nrow = 2, byrow = TRUE)
#' delta = c(0.5, 0.5)
#' # simulation
#' s = x = rep(NA, 500)
#' s[1] = sample(1:2, 1, prob = delta)
#' x[1] = rnorm(1, mu[s[1]], sigma[s[1]])
#' for(t in 2:500){
#'   s[t] = sample(1:2, 1, prob = Gamma[s[t-1],])
#'   x[t] = rnorm(1, mu[s[t]], sigma[s[t]])
#' }
#' 
#' ## negative log likelihood function
#' mllk = function(theta.star, x){
#'   # parameter transformations for unconstraint optimization
#'   Gamma = tpm(theta.star[1:2])
#'   delta = stationary(Gamma) # stationary HMM
#'   mu = theta.star[3:4]
#'   sigma = exp(theta.star[5:6])
#'   # calculate all state-dependent probabilities
#'   allprobs = matrix(1, length(x), 2)
#'   for(j in 1:2){ allprobs[,j] = dnorm(x, mu[j], sigma[j]) }
#'   # return negative for minimization
#'   -forward(delta, Gamma, allprobs)
#' }
#' 
#' ## fitting an HMM to the data
#' theta.star = c(-2,-2,0,5,log(2),log(3))
#' mod = stats::nlm(mllk, theta.star, x = x)
#'
forward = function(delta, Gamma, allprobs){
  forward_cpp_h(allprobs, delta, Gamma)
}
