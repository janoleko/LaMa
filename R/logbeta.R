logbeta = function(Gamma, allprobs){
  n = nrow(allprobs)
  N = ncol(allprobs)
  lbeta = matrix(NA, n, N)
  
  lbeta[n,] = rep(0,N)
  foo = rep(1/N, N)
  l = log(N)
  
  for(t in (n-1):1){
    foo = Gamma[,,t]%*%diag(allprobs[t+1,])%*%foo
    lbeta[t,] = l + log(foo)
    l = l + log(sum(foo))
    foo = foo/sum(foo)
  }
  return(lbeta)
}