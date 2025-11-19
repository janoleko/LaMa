# Computes generalised determinant

Computes generalised determinant

## Usage

``` r
gdeterminant(x, eps = NULL, log = TRUE)
```

## Arguments

- x:

  symmetric matrix

- eps:

  eigenvalues smaller than this will be treated as zero

- log:

  logical. If `TRUE`, the log-determinant is returned. If `FALSE`, the
  determinant is returned.

## Value

generalised log-determinant of `x`
