# AD-compatible minimum and maximum functions

These functions compute the parallel minimum/ maximum of two
vector-valued inputs and are compatible with automatic differentiation
using `RTMB`.

## Usage

``` r
min2(x, y)

max2(x, y)
```

## Arguments

- x:

  first vector

- y:

  second vector

## Value

`min2` returns the parallel minimum and `max2` the parallel maximum of
`x` and `y`

## Examples

``` r
x <- c(1, 4, 8, 2)
y <- c(2, 5, 3, 7)
min2(x, y)
#> [1] 1 4 3 2
max2(x, y)
#> [1] 2 5 8 7
```
