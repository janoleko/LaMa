logalpha = function(delta, Gamma, allprobs){
  n = nrow(allprobs)
  N = ncol(allprobs)
  lalpha = matrix(NA, n, N)
  
  foo = delta %*% diag(allprobs[1,])
  l = log(sum(foo))
  phi = foo/sum(foo)
  
  lalpha[1,] = l + log(foo)
  for(t in 2:n){
    foo = phi%*%Gamma[,,t-1]%*%diag(allprobs[t,])
    l = l + log(sum(foo))
    phi = foo/sum(foo)
    lalpha[t,] = l + log(foo)
  }
  return(lalpha)
}