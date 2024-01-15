
# {Lcpp}: Convenient forward algorithm in C++

This package contains convenient R-wrapper functions for the forward
algorithm used to fit hidden Markov models (HMMs) and state space models (SSMs) via direct
numerical maximum likelihood estimation. The algorithm calculates the log-likelihood recursively as a matrix product and uses a scaling strategy to avoid numerical underflow (see [Zucchini
et
al. 2016](https://www.taylorfrancis.com/books/mono/10.1201/b20790/hidden-markov-models-time-series-walter-zucchini-iain-macdonald-roland-langrock)).
Implementation in C++ offers 10-20 times faster evaluation times and thus substantially speeds up numerical maximization using e.g. `nlm()` or `optim()`.
Currently, the two main functions are `forward()` which can be used to fit homogeneous HMMs and `forward_g()` which can by used to fit inhomogeneous HMMs.
The functions are supposed to be included in the negative
log-likelihood function, once the parameter transformations and the allprobs matrix are calculated.

Further algorithm variations will be added as needed. Have fun!

## Installation

``` r
devtools::install_github("janoleko/Lcpp")
```
