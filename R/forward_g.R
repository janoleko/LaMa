#' General \href{https://www.taylorfrancis.com/books/mono/10.1201/b20790/hidden-markov-models-time-series-walter-zucchini-iain-macdonald-roland-langrock}{forward algorithm} with time-varying transition probability matrix
#'
#' @param delta Initial distribution of length N
#' @param Gamma Array of transition probability matrices of dimension c(N,N,n). \cr \cr
#' Here we use the definition \eqn{\Pr(S_t=j \mid S_{t-1}=i) = \gamma_{ij}^{(t)}}
#' such that the transition probabilities between time point \eqn{t-1} and \eqn{t} are an element of \eqn{\Gamma^{(t)}}.
#' Therefore, the first element of the array is not used in the likelihood calculation. \cr \cr
#' This function can also be used to fit continuous-time HMMs, where each array entry is the Markov semigroup \eqn{\Gamma(\Delta t) = \exp(Q \Delta t)} and \eqn{Q} is the generator.
#' @param allprobs Matrix of state-dependent probabilities/ density values of dimension c(n, N)
#'
#' @return Log-likelihood for given data and parameters
#' @export
#'
#' @examples
#' ## generating data from inhomogeneous 2-state HMM
#' mu = c(0, 6)
#' sigma = c(2, 4)
#' beta = matrix(c(-2,-2,0.5,-0.5),nrow=2)
#' delta = c(0.5, 0.5)
#' # simulation
#' n = 2000
#' s = x = rep(NA, n)
#' z = rnorm(n, 0, 2)
#' s[1] = sample(1:2, 1, prob = delta)
#' x[1] = stats::rnorm(1, mu[s[1]], sigma[s[1]])
#' for(t in 2:n){
#'   Gamma = diag(2)
#'   Gamma[!Gamma] = exp(beta[,1]+beta[,2]*z[t])
#'   Gamma = Gamma / rowSums(Gamma)
#'   s[t] = sample(1:2, 1, prob = Gamma[s[t-1],])
#'   x[t] = stats::rnorm(1, mu[s[t]], sigma[s[t]])
#' }
#' 
#' ## negative log likelihood function
#' mllk = function(theta.star, x, z){
#'   # parameter transformations for unconstraint optimization
#'   beta = matrix(theta.star[1:4], 2, 2)
#'   Eta = cbind(1, z)%*%t(beta)
#'   n = length(x)
#'   Gamma = array(dim=c(2,2,n))
#'   for(t in 2:n){
#'     G = diag(2)
#'     G[!G] = exp(Eta[t,])
#'     Gamma[,,t] = G / rowSums(G)
#'   }
#'   delta = c(plogis(theta.star[5]), 1-plogis(theta.star[5]))
#'   mu = theta.star[6:7]
#'   sigma = exp(theta.star[8:9])
#'   # calculate all state-dependent probabilities
#'   allprobs = matrix(1, length(x), 2)
#'   for(j in 1:2){ allprobs[,j] = stats::dnorm(x, mu[j], sigma[j]) }
#'   # return negative for minimization
#'   -Lcpp::forward_g(delta, Gamma, allprobs)
#' }
#' 
#' ## fitting an HMM to the data
#' theta.star = c(-2,-2,1,-1,0,0,5,log(2),log(3))
#' mod = stats::nlm(mllk, theta.star, x = x, z = z)
#'
forward_g = function(delta, Gamma, allprobs){
  forward_cpp_g(allprobs, delta, Gamma)
}