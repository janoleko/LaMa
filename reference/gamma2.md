# Reparametrised gamma distribution

Density, distribution function, quantile function and random generation
for the gamma distribution reparametrised in terms of mean and standard
deviation.

## Usage

``` r
dgamma2(x, mean = 1, sd = 1, log = FALSE)

pgamma2(q, mean = 1, sd = 1, lower.tail = TRUE, log.p = FALSE)

qgamma2(p, mean = 1, sd = 1, lower.tail = TRUE, log.p = FALSE)

rgamma2(n, mean = 1, sd = 1)
```

## Arguments

- x, q:

  vector of quantiles

- mean:

  mean parameter, must be positive scalar.

- sd:

  standard deviation parameter, must be positive scalar.

- log, log.p:

  logical; if `TRUE`, probabilities/ densities \\p\\ are returned as
  \\\log(p)\\.

- lower.tail:

  logical; if `TRUE`, probabilities are \\P\[X \<= x\]\\, otherwise,
  \\P\[X \> x\]\\.

- p:

  vector of probabilities

- n:

  number of observations. If `length(n) > 1`, the length is taken to be
  the number required.

## Value

`dgamma2` gives the density, `pgamma2` gives the distribution function,
`qgamma2` gives the quantile function, and `rgamma2` generates random
deviates.

## Details

This implementation allows for automatic differentiation with `RTMB`.

## Examples

``` r
x = rgamma2(1)
d = dgamma2(x)
p = pgamma2(x)
q = qgamma2(p)
```
