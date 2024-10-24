
# LaMa <img src="man/figures/Logo_LaMa.png" align="right" height=150>

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/LaMa)](https://CRAN.R-project.org/package=LaMa)
[![metacran
downloads](https://cranlogs.r-pkg.org/badges/last-month/LaMa)](https://cran.r-project.org/package=LaMa)
[![R-CMD-check](https://github.com/janoleko/LaMa/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/janoleko/LaMa/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

A variety of latent Markov models, including **hidden Markov models**
(HMMs), **hidden semi-Markov models** (HSMMs), **state space models**
(SSMs) and **continuous-time** variants can be formulated and estimated
within the same framework via directly maximizing the (approximate)
likelihood using the so-called **forward algorithm** (for details see
<a href="https://www.taylorfrancis.com/books/mono/10.1201/b20790/hidden-markov-models-time-series-walter-zucchini-iain-macdonald-roland-langrock" target="_blank">Zucchini
et al. 2016</a>). Applied researchers often need custom models that
standard software does not easily support. Writing tailored `R` code
offers flexibility but suffers from slow estimation speeds. This `R`
package solves these issues by providing easy-to-use functions (written
in C++ for speed) for common tasks like the forward algorithm. These
functions can be combined into custom models, offering up to 10-20 times
faster estimation via standard numerical optimizers like `nlminb()` or
`optim()`. The development version now also allows for automatic
differentiation with the `RTMB` package which drastically increases
speed and accuracy.

The most important implementations of the forward algorithm are:

- `forward()` for models with **homogeneous** transition probabilities,
  and
- `forward_g()` for general (pre-calculated) **inhomogeneous**
  transition probabilities (including **continuous-time** HMMs and
  points processes)

The functions are built to be included in the **negative log-likelihood
function**, after parameters have been transformed and the `allprobs`
matrix (containing all state-dependent probabilities) has been
calculated.

To serve as a powerful toolbox, this package also includes many
auxiliary functions like

- the `tpm` family with for calculating transition probability matrices,
  <!-- + `tpm()` for calculating a homogeneous transition probability matrix via the multinomial logistic link,  -->
  <!-- + `tpm_g()` for calculating general inhomogeneous transition probabilty matrices,  -->
  <!-- + `tpm_p()` for calculating transition matrices of periodically inhomogeneous HMMs, -->
  <!-- + `tpm_cont()` for calculating the transition probabilites of a continuous-time Markov chain, -->
  <!-- + `tpm_emb()` for calculating the embedded transition matrix of an HSMM, -->

- the `stationary` family to compute stationary and periodically
  stationary distributions

- as well as the `stateprobs` and `viterbi` family for local and global
  decoding.

## Installation

You can install the released package version from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("LaMa")
```

or the development version from Github:

``` r
# install.packages("devtools")
devtools::install_github("janoleko/LaMa")
```

<!-- (To install from Github, you need a functional <a href="https://teuder.github.io/rcpp4everyone_en/020_install.html" target="_blank">C++ compiler</a>.) -->

## Package documentation

To aid in building fully custom likelihood functions, this package
contains several vignettes that demonstrate how to simulate data from
and estimate a wide range of models using the functions included in this
package:

- [Introduction to
  LaMa](https://janoleko.github.io/LaMa/articles/Intro_to_LaMa.html)
- [Inhomogeneous HMMs with covariate
  effects](https://janoleko.github.io/LaMa/articles/Inhomogeneous_HMMs.html)
- [Longitudinal
  data](https://janoleko.github.io/LaMa/articles/Longitudinal_data.html)
- [Periodic
  HMMs](https://janoleko.github.io/LaMa/articles/Periodic_HMM.html)
- [State space
  models](https://janoleko.github.io/LaMa/articles/State_space_models.html)
- [Continuous-time
  HMMs](https://janoleko.github.io/LaMa/articles/Continuous_time_HMMs.html)
- [Hidden semi-Markov
  models](https://janoleko.github.io/LaMa/articles/HSMMs.html)
- [Markov-modulated (marked) Poisson
  processes](https://janoleko.github.io/LaMa/articles/MMMPPs.html)
- [LaMa and
  RTMB](https://janoleko.github.io/LaMa/articles/LaMa_and_RTMB.html)

<!-- ## Citation -->
<!-- When using LaMa, please cite the package as follows: -->
<!-- ```{r citation} -->
<!-- citation(package = "LaMa") -->
<!-- ``` -->

## Introductory example: Homogeneous HMM

We analyze the `elephant` data set contained in the package using a
simple 2-state HMM with state-dependent gamma distributions. We fit the
model to the hourly step lengths.

``` r
library(LaMa)
#> Loading required package: RTMB

head(elephant, 3)
#>   tod      step     angle state
#> 1   9 0.3252437        NA     1
#> 2  10 0.2458265  2.234562     1
#> 3  11 0.2173252 -2.262418     1
```

We start by defining the negative log-likelihood function. This is made
really convenient by the functions `tpm()` which computes the transition
probability matrix via the multinomial logit link, `stationary()` which
computes the stationary distribution of the Markov chain and `forward()`
which calculates the log-likelihood via the forward algorithm.

``` r
nll = function(par, step){
  # parameter transformations for unconstraint optimization
  Gamma = tpm(par[1:2]) # multinomial logit link
  delta = stationary(Gamma) # stationary HMM
  mu = exp(par[3:4])
  sigma = exp(par[5:6])
  # calculate all state-dependent probabilities
  allprobs = matrix(1, length(step), 2)
  ind = which(!is.na(step))
  for(j in 1:2) allprobs[ind,j] = dgamma2(step[ind], mu[j], sigma[j])
  # simple forward algorithm to calculate log-likelihood
  -forward(delta, Gamma, allprobs)
}
```

To fit the model, we define the intial parameter vector and numerically
optimize the above function using `nlm()`:

``` r
par = c(-2,-2,             # initial tpm params (logit-scale)
        log(c(0.3, 1)),    # initial means for step length (log-transformed)
        log(c(0.2, 0.7)))  # initial sds for step length (log-transformed)

# initial transformed parameters: not chosen too well
system.time(
  mod <- nlm(nll, par, step = elephant$step)
)
#>    user  system elapsed 
#>   0.456   0.012   0.468
```

Really fast for 10.000 data points!

After tranforming the unconstrained parameters to working parameters
using `tpm()` and `stationary()`, we can visualize the results:

``` r
# transform parameters to working
Gamma = tpm(mod$estimate[1:2])
delta = stationary(Gamma) # stationary HMM
mu = exp(mod$estimate[3:4])
sigma = exp(mod$estimate[5:6])

hist(elephant$step, prob = TRUE, bor = "white", breaks = 40, main = "", xlab = "step length")
curve(delta[1] * dgamma2(x, mu[1], sigma[1]), add = TRUE, lwd = 2, col = "orange", n=500)
curve(delta[2] * dgamma2(x, mu[2], sigma[2]), add = TRUE, lwd = 2, col = "deepskyblue", n=500)
legend("topright", col = c("orange", "deepskyblue"), lwd = 2, bty = "n", legend = c("state 1", "state 2"))
```

<img src="man/figures/README-visualization-1.png" width="75%" style="display: block; margin: auto;" />
