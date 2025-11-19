# Zero-inflated density constructer

Constructs a zero-inflated density function from a given probability
density function

## Usage

``` r
zero_inflate(dist, discrete = NULL)
```

## Arguments

- dist:

  either a probability density function or a probability mass function

- discrete:

  logical; if `TRUE`, the density for `x = 0` will be
  `zeroprob + (1-zeroprob) * dist(0, ...)`. Otherwise it will just be
  `zeroprob`. In standard cases, this will be determined automatically.
  For non-standard cases, set this to `TRUE` or `FALSE` depending on the
  type of `dist`. See details.

## Value

zero-inflated density function with first argument `x`, second argument
`zeroprob`, and additional arguments `...` that will be passed to
`dist`.

## Details

The definition of zero-inflation is different for discrete and
continuous distributions. For discrete distributions with p.m.f. \\f\\
and zero-inflation probability \\p\\, we have \$\$\Pr(X = 0) = p + (1 -
p) \cdot f(0),\$\$ and \$\$\Pr(X = x) = (1 - p) \cdot f(x), \quad x \>
0.\$\$

For continuous distributions with p.d.f. \\f\\, we have
\$\$f\_{\text{zinfl}}(x) = p \cdot \delta_0(x) + (1 - p) \cdot f(x),\$\$
where \\\delta_0\\ is the Dirac delta function at zero.

## Examples

``` r
dzinorm <- zero_inflate(dnorm)
dzinorm(c(NA, 0, 2), 0.5, mean = 1, sd = 1)
#> [1]        NA 0.5000000 0.1209854

zipois <- zero_inflate(dpois)
zipois(c(NA, 0, 1), 0.5, 1)
#> [1]        NA 0.6839397 0.1839397
```
