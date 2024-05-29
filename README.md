
# {LaMa}: Latent Markov model likelihood evaluation in C++ <img src="man/figures/Logo_LaMa.png" align="right" height=170>

A plethora of latent Markov models, including **hidden Markov models**
(HMMs), **hidden semi-Markov models** (HSMMs), **state space models**
(SSMs) as well as **continuous-time HMMs**, **continuous-time SSMs**,
and **Markov-modulated marked Poisson processes** (MMMPPs) can be
formulated and estimated within the same framework via directly
maximizing the (approximate) likelihood using the so-called **forward
algorithm** (for details see
<a href="https://teuder.github.io/rcpp4everyone_en/020_install.html" target="_blank">Zucchini
et al. 2016</a>). While these models are immensely popular for
conducting inference on time series driven by latent processes, due to
their great flexibility, researchers using these models in applied work
often need to build highly customized models for which standard software
implementation is lacking, or the construction of such models in said
software is as complicated as writing fully tailored **R** code. The
latter provides great flexibility and control, but suffers from slow
estimation speeds that make custom solutions inconvenient. This **R**
package addresses the above issues in two ways. Standard blocks of code
common to all these model classes, most importantly the forward
algorithm, are implemented as simple-to-use functions. These can be
added like Lego blocks to an otherwise fully custom likelihood function,
making building fully custom models much easier. Moreover, under the
hood, these functions are written in **C++**, allowing for 10-20 times
faster evaluation time, and thus drastically speeding up estimation by
numerical optimizers like `nlm()` or `optim()`.

Current implementations of the forward algorithm are:

- `forward()` for models with **homogeneous** transition probabilities,
- `forward_g()` for general (pre-calculated) **inhomogeneous**
  transition probabilities (including **continuous-time** HMMs and
  points processes), and
- `forward_s()` for fitting **HSMMs**.

The functions are built to be included in the **negative log-likelihood
function**, after parameters have been transformed and the *allprobs*
matrix (containing all state-dependent probabilities) has been
calculated.

To serve as a powerful toolbox, this package also includes many
auxiliary functions like

- The `tpm` family with

  - `tpm()` for calculating a homogeneous transition probability matrix
    via the multinomial logistic link,
  - `tpm_g()` for calculating general inhomogeneous transition
    probabilty matrices,
  - `tpm_p()` for calculating transition matrices of periodically
    inhomogeneous HMMs,
  - `tpm_cont()` for calculating the transition probabilites of a
    continuous-time Markov chain,
  - `tpm_hsmm()` for calculating the transition matrix of an
    HSMM-approximating HMM,

- the `stationary` family to compute stationary and periodically
  stationary distributions,

- functions of the `stateprobs` family for local decoding and of the
  `viterbi` family for global decoding,

- and `trigBasisExp()` for efficient computation of a trigonometric
  basis expansion.

Further functionalities will be added as needed. Have fun!

## Package documentation

To aid in building fully custom likelihood functions, this package also
contains several vignettes that show how to simulate data from and
estimate a wide range of models:

- [Introduction to
  LaMa](https://janoleko.github.io/files/vignettes/LaMa/Intro_to_LaMa.pdf)
- [Inhomogeneous
  HMMs](https://janoleko.github.io/files/vignettes/LaMa/Inhomogeneous_HMMs.pdf)
- [Longitudinal
  data](https://janoleko.github.io/files/vignettes/LaMa/Longitudinal_data.pdf)
- [Periodic
  HMMs](https://janoleko.github.io/files/vignettes/LaMa/Periodic_HMM.pdf)
- [State space
  models](https://janoleko.github.io/files/vignettes/LaMa/State_space_models.pdf)
- [Continuous-time
  HMMs](https://janoleko.github.io/files/vignettes/LaMa/Continuous_time_HMMs.pdf)
- [Hidden semi-Markov
  models](https://janoleko.github.io/files/vignettes/LaMa/HSMMs.pdf)
- [Markov-modulated (marked) Poisson
  processes](https://janoleko.github.io/files/vignettes/LaMa/MMMPPs.pdf)

## Installation

To install and use the package, you need to have a functional C++
compiler. For details click
<a href="https://teuder.github.io/rcpp4everyone_en/020_install.html" target="_blank">here</a>.
Then you can use:

``` r
# install.packages("devtools")
devtools::install_github("janoleko/LaMa", build_vignettes = TRUE)
```

Feel free to use

``` r
browseVignettes("LaMa")
```

or

``` r
help(package = "LaMa")
```

for detailed examples on how to use the package.

## Example: Homogeneous HMM

#### Loading the package

``` r
library(LaMa)
```

#### Generating data from a 2-state HMM

Here we can use `stationary()` to compute the stationary distribution.

``` r
# parameters
mu = c(0, 6)
sigma = c(2, 4)
Gamma = matrix(c(0.95, 0.05, 0.15, 0.85), nrow = 2, byrow = TRUE)
delta = stationary(Gamma) # stationary HMM

# simulation
n = 10000 # rather large
set.seed(123)
s = x = rep(NA, n)
s[1] = sample(1:2, 1, prob = delta)
x[1] = rnorm(1, mu[s[1]], sigma[s[1]])
for(t in 2:n){
  s[t] = sample(1:2, 1, prob = Gamma[s[t-1],])
  x[t] = rnorm(1, mu[s[t]], sigma[s[t]])
}

plot(x[1:200], bty = "n", pch = 20, ylab = "x", 
     col = c("orange","deepskyblue")[s[1:200]])
```

<img src="man/figures/README-data-1.png" width="75%" style="display: block; margin: auto;" />

#### Writing the negative log-likelihood function

Here, we build the transition probability matrix using the `tpm()`
function, compute the stationary distribution using `stationary()` and
calculate the log-likelihood using `forward()` in the last line.

``` r
mllk = function(theta.star, x){
  # parameter transformations for unconstraint optimization
  Gamma = tpm(theta.star[1:2])
  delta = stationary(Gamma) # stationary HMM
  mu = theta.star[3:4]
  sigma = exp(theta.star[5:6])
  # calculate all state-dependent probabilities
  allprobs = matrix(1, length(x), 2)
  for(j in 1:2){ allprobs[,j] = stats::dnorm(x, mu[j], sigma[j]) }
  # return negative for minimization
  -forward(delta, Gamma, allprobs)
}
```

#### Fitting an HMM to the data

``` r
theta.star = c(-1,-1,1,4,log(1),log(3)) 
# initial transformed parameters: not chosen too well
s = Sys.time()
mod = nlm(mllk, theta.star, x = x)
Sys.time()-s
#> Time difference of 0.1056249 secs
```

Really fast for 10.000 data points!

#### Visualizing results

Again, we use `tpm()` and `stationary()` to tranform the unconstraint
parameters to working parameters.

``` r
# transform parameters to working
Gamma = tpm(mod$estimate[1:2])
delta = stationary(Gamma) # stationary HMM
mu = mod$estimate[3:4]
sigma = exp(mod$estimate[5:6])

hist(x, prob = TRUE, bor = "white", breaks = 40, main = "")
curve(delta[1]*dnorm(x, mu[1], sigma[1]), add = TRUE, lwd = 2, col = "orange", n=500)
curve(delta[2]*dnorm(x, mu[2], sigma[2]), add = TRUE, lwd = 2, col = "deepskyblue", n=500)
curve(delta[1]*dnorm(x, mu[1], sigma[1])+delta[2]*dnorm(x, mu[2], sigma[2]),
      add = TRUE, lwd = 2, lty = "dashed", n=500)
legend("topright", col = c("orange", "deepskyblue", "black"), lwd = 2, bty = "n",
       lty = c(1,1,2), legend = c("state 1", "state 2", "marginal"))
```

<img src="man/figures/README-visualization-1.png" width="75%" style="display: block; margin: auto;" />
