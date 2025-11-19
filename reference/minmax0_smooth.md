# Smooth approximations to max(x, 0) and min(x, 0)

Smooth approximations to max(x, 0) and min(x, 0)

## Usage

``` r
max0_smooth(x, rho = 20)

min0_smooth(x, rho = 20)
```

## Arguments

- x:

  a vector of values

- rho:

  smoothing parameter, larger values lead to closer approximation

## Value

the approximate maximum or minimum of x and 0

## Examples

``` r
x <- seq(-1, 1, by = 0.1)
min0_smooth(x)
#>  [1] -1.000000e+00 -9.000000e-01 -8.000000e-01 -7.000000e-01 -6.000003e-01
#>  [6] -5.000023e-01 -4.000168e-01 -3.001238e-01 -2.009075e-01 -1.063464e-01
#> [11] -3.465736e-02 -6.346401e-03 -9.074964e-04 -1.237843e-04 -1.677032e-05
#> [16] -2.269945e-06 -3.072097e-07 -4.157642e-08 -5.626758e-09 -7.614990e-10
#> [21] -1.030577e-10
max0_smooth(x)
#>  [1] 1.030577e-10 7.614990e-10 5.626758e-09 4.157642e-08 3.072097e-07
#>  [6] 2.269945e-06 1.677032e-05 1.237843e-04 9.074964e-04 6.346401e-03
#> [11] 3.465736e-02 1.063464e-01 2.009075e-01 3.001238e-01 4.000168e-01
#> [16] 5.000023e-01 6.000003e-01 7.000000e-01 8.000000e-01 9.000000e-01
#> [21] 1.000000e+00
```
