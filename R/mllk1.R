#' Negative log Likelihood for the fly data RE model
#'
#' @param theta.star unconstraint working parameter vector
#' @param X data
#' @param L number of time points in a cycle
#' @param K degree of trigonometric link
#' @param M numer of intervals for numerical approximation of the integral
#' @param low lower integral boundaries
#' @param up upper integral boundaries
#'
#' @return returns the negative log likelihood of the data given the paramter vector
#' @export
#'
#' @examples
mllk_random = function(theta.star, X, L = 48, K, M, low, up){
  N = 2
  coef = array(NA, dim = c(N*(N-1), 1+2*K, 2))
  coef[,,1] = matrix(theta.star[1:((1+2*K)*N*(N-1))], nrow = N*(N-1), ncol = (1+2*K))
  coef[,,2] = matrix(theta.star[((1+2*K)*(N-1)*N) + 1:((1+2*K)*N*(N-1))], nrow = N*(N-1), ncol = (1+2*K))
  mu = exp(theta.star[(2*(1+2*K)*(N-1)*N) + 1:N])
  sigma = exp(theta.star[(2*(1+2*K)*(N-1)*N) + N + 1:N])
  size = exp(theta.star[(2*(1+2*K)*(N-1)*N) + 2*N + 1:N])
  Gamma = array(NA, dim = c(N,N,L,2))
  for (k in 1:L){
    G = diag(N)
    G[!G] = exp(iHSMM::pv(coef[,,1], time = (k-1)/2, degree = K, L = L))
    G = G/rowSums(G)
    Gamma[,,k,1] = G
    G = diag(N)
    G[!G] = exp(iHSMM::pv(coef[,,2], time = (k-1)/2, degree = K, L = L))
    G = G/rowSums(G)
    Gamma[,,k,2] = G }
  GammaT=Gamma[,,X$tod[1],1]
  for(k in 42+1:(L-1)){
    k=ifelse(k%%L==0,L,k%%L)
    GammaT = GammaT%*%Gamma[,,k,1] }
  delta = solve(t(diag(nrow(GammaT))-GammaT+1), rep(1,nrow(GammaT)))
  a1 = seq(from=low[1], to=up[1], length.out = M) # Intervals for first RE
  a2 = seq(from=low[2], to=up[2], length.out = M) # Intervals for second RE
  a.m1 = (a1[-1] + a1[-M])*0.5 # midpoints interval first RE
  a.m2 = (a2[-1] + a2[-M])*0.5 # midpoints interval second RE 
  am = cbind(a.m1, a.m2)
  stan.fac = (stats::pgamma(max(a1), shape=mu[1]^2/sigma[1]^2,scale=sigma[1]^2/mu[1])-stats::pgamma(min(a1), shape=mu[1]^2/sigma[1]^2,scale=sigma[1]^2/mu[1]))*
    (stats::pgamma(max(a2), shape=mu[2]^2/sigma[2]^2,scale=sigma[2]^2/mu[2])-stats::pgamma(min(a2), shape=mu[2]^2/sigma[2]^2,scale=sigma[2]^2/mu[2]))
  probs = numeric((M-1)^2)
  for(j1 in 1:(M-1)){
    for(j2 in 1:(M-1)){
      r = (j1-1)*(M-1)+1+(j2-1)
      probs[r] = (stats::pgamma(a1[j1+1], shape=mu[1]^2/sigma[1]^2,scale=sigma[1]^2/mu[1])-stats::pgamma(a1[j1], shape=mu[1]^2/sigma[1]^2,scale=sigma[1]^2/mu[1]))*
        (stats::pgamma(a2[j2+1], shape=mu[2]^2/sigma[2]^2,scale=sigma[2]^2/mu[2])-stats::pgamma(a2[j2], shape=mu[2]^2/sigma[2]^2,scale=sigma[2]^2/mu[2]))
    }}
  probs = probs/stan.fac
  IDs = unique(X$ID)
  nbAnimals = length(IDs)
  l = 0
  for(k in 1:nbAnimals){
    X_k = X[which(X$ID == IDs[k]),]
    startDD = which(X_k$condition == "DD")[1]
    nObs = nrow(X_k)
    ind = which(!is.na(X_k$activity))
    one.animal = numeric((M-1)^2)
    for(j1 in 1:(M-1)){
      for(j2 in 1:(M-1)){
        allprobs = matrix(1, nObs, N)
        for (j in 1:N){ 
          eta = c(am[j1,1], am[j2,2])
          allprobs[ind,j] = stats::dnbinom(X_k$activity[ind], size = size[j], mu = eta[j]) }
        llk = forward_cpp(allprobs, delta, Gamma[,,,1], Gamma[,,,2], startDD, X_k$tod-1)
        r = (j1-1)*(M-1)+1+(j2-1)
        one.animal[r] = llk
      }}   
    ma = max(one.animal)
    logL = numeric(nObs)
    notNull = which(ma - one.animal < 700)
    logL[notNull] = exp(one.animal[notNull] - ma)
    l = l + (log(sum(logL*probs)) + ma) 
  }
  return(-l)
}