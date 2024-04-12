#' General \href{https://www.taylorfrancis.com/books/mono/10.1201/b20790/hidden-markov-models-time-series-walter-zucchini-iain-macdonald-roland-langrock}{forward algorithm} with time-varying transition probability matrix
#'
#' @param delta Initial distribution of length N, or matrix of dimension c(N,k) for k independent tracks, if trackInd is provided.
#' @param Gamma Array of transition probability matrices of dimension c(N,N,n-1), as in a time series of length n, there are only n-1 transitions. 
#' If you provide an array of dimension c(N,N,n), the first slice will be ignored. \cr
#' 
#' If the elements of \eqn{\Gamma^{(t)}} depend on covariate values at t or covariates t+1 is your choice in the calculation of the array, prior to using this function.
#' When conducting the calculation by using tpm_g(), the choice comes down to including the covariate matrix Z[-1,] oder Z[-n,]. \cr
#' 
#' This function can also be used to fit continuous-time HMMs, where each array entry is the Markov semigroup \eqn{\Gamma(\Delta t) = \exp(Q \Delta t)} and \eqn{Q} is the generator of the continuous-time Markov chain.
#' 
#' @param allprobs Matrix of state-dependent probabilities/ density values of dimension c(n, N)
#' @param trackInd Optional vector of length k containing the indices that correspond to the beginning of a separate track. If provided, the total log-likelihood will be the sum of each track's likelihood contribution.
#' In this case, Gamma needs to be an array of dimension c(N,N,n), matching the number of rows of allprobs. For each track, the transition matrix at the beginning of the track will be ignored (as there is no transition between tracks).
#' Furthermore, instead of a single vector delta corresponding to the initial distribution, a delta matrix of initial distributions, of dimension c(N,k), can be provided, such that each track starts with it's own initial distribution.
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
#' x[1] = rnorm(1, mu[s[1]], sigma[s[1]])
#' for(t in 2:n){
#'   Gamma = diag(2)
#'   Gamma[!Gamma] = exp(beta[,1]+beta[,2]*z[t])
#'   Gamma = Gamma / rowSums(Gamma)
#'   s[t] = sample(1:2, 1, prob = Gamma[s[t-1],])
#'   x[t] = rnorm(1, mu[s[t]], sigma[s[t]])
#' }
#' 
#' ## negative log likelihood function
#' mllk = function(theta.star, x, z){
#'   # parameter transformations for unconstraint optimization
#'   beta = matrix(theta.star[1:4], 2, 2)
#'   Gamma = tpm_g(Z = z, beta = beta)
#'   delta = c(plogis(theta.star[5]), 1-plogis(theta.star[5]))
#'   mu = theta.star[6:7]
#'   sigma = exp(theta.star[8:9])
#'   # calculate all state-dependent probabilities
#'   allprobs = matrix(1, length(x), 2)
#'   for(j in 1:2){ allprobs[,j] = dnorm(x, mu[j], sigma[j]) }
#'   # return negative for minimization
#'   -forward_g(delta, Gamma, allprobs)
#' }
#' 
#' ## fitting an HMM to the data
#' theta.star = c(-2,-2,1,-1,0,0,5,log(2),log(3))
#' mod = nlm(mllk, theta.star, x = x, z = z)
#'
forward_g = function(delta, Gamma, allprobs, trackInd = NULL){
  n = nrow(allprobs)
  
  if(is.null(trackInd)){
    if(dim(Gamma)[3]==n){
      Gamma = Gamma[,,-1]
    }
    l = forward_cpp_g(allprobs, delta, Gamma)
  } else{
    k = length(trackInd)
    
    if(dim(Gamma)[3]!=n) stop("Gamma needs to be an array of dimension c(N,N,n), matching the number of rows of allprobs.")
    
    if(is.vector(delta)){
      delta = matrix(delta, nrow = k, ncol = length(delta))
    }
    
    l = forward_cpp_g_tracks(allprobs, delta, Gamma, trackInd)
  }
  
  return(l)
}
