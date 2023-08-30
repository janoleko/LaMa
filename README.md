
# {Lcpp}: Convenient forward algorithm in C++

This package contains convenient R-wrapper functions of the forward
algorithm used to estimate hidden Markov models (HMMs) via direct
numerical maximum likelihood estimation. The algorithm calculates the
likelihood recursively as a matrix product (see [Zucchini et
al. 2016](https://www.taylorfrancis.com/books/mono/10.1201/b20790/hidden-markov-models-time-series-walter-zucchini-iain-macdonald-roland-langrock)).
Thus, implemntation in C++ offers 10-20 times faster evaluation times.
Several different functions are contained including homogeneous HMMs and
a general implementation with a pre-calculated array of transition
probability matrices (t.p.m.s). These functions can then be included in
the negative log-likelihood function, once the initial distribution
(delta), the allprobs matrix and the t.p.m.(s) are calculated.

Various different implementations will be added as needed. Have fun!

## Installation

``` r
devtools::install_github("janoleko/Lcpp")
```
