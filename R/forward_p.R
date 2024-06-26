#' \href{https://www.taylorfrancis.com/books/mono/10.1201/b20790/hidden-markov-models-time-series-walter-zucchini-iain-macdonald-roland-langrock}{Forward algorithm} with (only) periodically varying transition probability matrix
#'
#' When the transition probability matrix only varies periodically (e.g. as a function of time of day), there are only \eqn{L} unique matrices if \eqn{L} is the period length (e.g. \eqn{L=24} for hourly data and time-of-day variation).
#' Thus, it is much more efficient to only calculate these \eqn{L} matrices and index them by a time variable (e.g. time of day or day of year) instead of calculating such a matrix for each index in the data set (which would be redundant).
#' This function allows for that, by only expecting a transition probability matrix for each time point in a period, and an integer valued (\eqn{1, \dots, L}) time variable that maps the data index to the according time.
#'
#' @param delta Initial or periodically stationary distribution of length N
#' @param Gamma Array of transition probability matrices of dimension c(N,N,L). \cr \cr
#' Here we use the definition \eqn{\Pr(S_t=j \mid S_{t-1}=i) = \gamma_{ij}^{(t)}}
#' such that the transition probabilities between time point \eqn{t-1} and \eqn{t} are an element of \eqn{\Gamma^{(t)}}.
#' @param allprobs Matrix of state-dependent probabilities/ density values of dimension c(n, N)
#' @param tod (Integer valued) time variable in 1, ..., L, mapping the data index to a generalized time of day (length n).
#' For half-hourly data L = 48. It could, however, also be day of year for daily data and L = 365.
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
#' x[1] = rnorm(1, mu[s[1]], sigma[s[1]])
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
#'   x[t] = rnorm(1, mu[s[t]], sigma[s[t]])
#' }
#' # we can also use function from LaMa to make building periodic tpms much easier
#' Gamma = tpm_p(1:24, 24, beta, degree = 1)
#' 
#' ## negative log likelihood function
#' mllk = function(theta.star, x, tod){
#'   # parameter transformations for unconstraint optimization
#'   beta = matrix(theta.star[1:6], 2, 3)
#'   Gamma = tpm_p(tod=tod, L=24, beta=beta)
#'   delta = stationary_p(Gamma, t=tod[1])
#'   mu = theta.star[8:9]
#'   sigma = exp(theta.star[10:11])
#'   # calculate all state-dependent probabilities
#'   allprobs = matrix(1, length(x), 2)
#'   for(j in 1:2){ allprobs[,j] = dnorm(x, mu[j], sigma[j]) }
#'   # return negative for minimization
#'   -forward_p(delta, Gamma, allprobs, tod)
#' }
#' 
#' ## fitting an HMM to the data
#' theta.star = c(-2,-2,1,-1,1,-1,0,0,5,log(2),log(3))
#' mod = nlm(mllk, theta.star, x = x, tod = tod)
#'
forward_p = function(delta, Gamma, allprobs, tod){
  if(min(tod)==1){
    tod = tod-1
  } 
  forward_cpp_p(allprobs, delta, Gamma, tod)
}