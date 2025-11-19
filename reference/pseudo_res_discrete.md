# Calculate pseudo-residuals for discrete-valued observations

For HMMs, pseudo-residuals are used to assess the goodness-of-fit of the
model. These are based on the cumulative distribution function (CDF)
\$\$F\_{X_t}(x_t) = F(x_t \mid x_1, \dots, x\_{t-1}, x\_{t+1}, \dots,
x_T)\$\$ and can be used to quantify whether an observation is extreme
relative to its model-implied distribution.

This function calculates such residuals for **discrete-valued**
observations, based on the local state probabilities obtained by
[`stateprobs`](https://janoleko.github.io/reference/stateprobs.md) or
[`stateprobs_g`](https://janoleko.github.io/reference/stateprobs_g.md)
and the respective parametric family.

## Usage

``` r
pseudo_res_discrete(
  obs,
  dist,
  par,
  stateprobs,
  normal = TRUE,
  randomise = TRUE,
  seed = NULL
)
```

## Arguments

- obs:

  vector of discrete-valued observations (of length n)

- dist:

  character string specifying which parametric CDF to use (e.g.,
  `"norm"` for normal or `"pois"` for Poisson) or CDF function to
  evaluate directly.

- par:

  named parameter list for the parametric CDF

  Names need to correspond to the parameter names in the specified
  distribution (e.g. `list(mean = c(1,2), sd = c(1,1))` for a normal
  distribution and 2 states). This argument is as flexible as the
  parametric distribution allows. For example you can have a matrix of
  parameters with one row for each observation and one column for each
  state.

- stateprobs:

  matrix of local state probabilities for each observation (of dimension
  c(n,N), where N is the number of states)

- normal:

  logical, if `TRUE`, returns Gaussian pseudo residuals

  These will be approximately standard normally distributed if the model
  is correct.

- randomise:

  logical, if `TRUE`, return randomised pseudo residuals. Recommended
  for discrete observations.

- seed:

  integer, seed for random number generation

## Value

vector of pseudo residuals

## Details

For discrete observations, calculating pseudo residuals is slightly more
involved, as the CDF is a step function. Therefore, one can calculate
the lower and upper CDF values for each observation. By default, this
function does exactly that and then randomly samples the interval in
between to give approximately Gaussian psuedo-residuals. If `randomise`
is set to `FALSE`, the lower, upper and mean pseudo-residuasl are
returned.

## Examples

``` r
obs = rpois(100, lambda = 1)
stateprobs = matrix(0.5, nrow = 100, ncol = 2)
par = list(lambda = c(1,2))
pres = pseudo_res_discrete(obs, "pois", par, stateprobs)
#> Randomised between lower and upper
```
