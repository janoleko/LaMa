
# {Lcpp}: Convenient forward algorithm in C++

This package contains convenient R-wrapper functions of the forward
algorithm used to estimate hidden Markov models (HMMs) via direct
numerical maximum likelihood estimation. The algorithm calculates the
negative log-likelihood recursively as a matrix product (see [Zucchini
et
al. 2016](https://www.taylorfrancis.com/books/mono/10.1201/b20790/hidden-markov-models-time-series-walter-zucchini-iain-macdonald-roland-langrock)).
Thus, implementation in C++ offers 10-20 times faster evaluation times.
Several versions are contained, including homogeneous HMMs and a general
implementation with a pre-calculated array of transition probability
matrices (t.p.m.s). These functions can then be included in the negative
log-likelihood function, once the initial distribution (delta), the
allprobs matrix and the t.p.m.(s) are calculated.

Various different implementations will be added as needed. Have fun!

## Installation

``` r
devtools::install_github("janoleko/Lcpp")
```

## Example

#### Generating some data from a 2-state HMM

``` r
# parameters
mu = c(0, 6)
sigma = c(2, 4)
Gamma = matrix(c(0.5, 0.05, 0.15, 0.85), nrow = 2, byrow = TRUE)
delta = c(0.5, 0.5)

# simulation
s = x = rep(NA, 500)
s[1] = sample(1:2, 1, prob = delta)
x[1] = stats::rnorm(1, mu[s[1]], sigma[s[1]])
for(t in 2:500){
  s[t] = sample(1:2, 1, prob = Gamma[s[t-1],])
  x[t] = stats::rnorm(1, mu[s[t]], sigma[s[t]])
}

plot(x, bty = "n", pch = 20, col = c("orange", "deepskyblue")[s])
```

<img src="man/figures/README-unnamed-chunk-3-1.png" width="75%" style="display: block; margin: auto;" />
\#### Writing the negative log-likelihood function

``` r
mllk = function(theta.star, x){
  # parameter transformations for unconstraint optimization
  Gamma = diag(2)
  Gamma[!Gamma] = exp(theta.star[1:2])
  Gamma = Gamma / rowSums(Gamma)
  delta = solve(t(diag(2)-Gamma+1), rep(1,2)) # stationary HMM
  mu = theta.star[3:4]
  sigma = exp(theta.star[5:6])
  # calculate all state-dependent probabilities
  allprobs = matrix(1, length(x), 2)
  for(j in 1:2){ allprobs[,j] = stats::dnorm(x, mu[j], sigma[j]) }
  # return negative for minimization
  -Lcpp::forward(delta, Gamma, allprobs)
}
```

#### Fit an HMM to the data

``` r
theta.star = c(-2,-2,0,5,log(2),log(3))
mod = stats::nlm(mllk, theta.star, x = x)
```

#### Visulize some results

``` r
# transform parameters
Gamma = diag(2)
Gamma[!Gamma] = exp(mod$estimate[1:2])
Gamma = Gamma / rowSums(Gamma)
delta = solve(t(diag(2)-Gamma+1), rep(1,2)) # stationary
mu = mod$estimate[3:4]
sigma = exp(mod$estimate[5:6])

hist(x, prob = T, bor = "white", breaks = 40, main = "Marginal distributions")
curve(delta[1]*dnorm(x, mu[1], sigma[1]), add = T, lwd = 2, col = "orange")
curve(delta[2]*dnorm(x, mu[2], sigma[2]), add = T, lwd = 2, col = "deepskyblue")
```

<img src="man/figures/README-unnamed-chunk-6-1.png" width="75%" style="display: block; margin: auto;" />
