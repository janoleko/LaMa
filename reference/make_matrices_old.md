# Build the design and the penalty matrix for models involving penalised splines based on a formula and a data set

Build the design and the penalty matrix for models involving penalised
splines based on a formula and a data set

## Usage

``` r
make_matrices_old(formula, data, knots = NULL)
```

## Arguments

- formula:

  right side of a formula as used in `mgcv`

- data:

  data frame containing the variables in the formula

- knots:

  optional list containing user specified knot values to be used for
  basis construction

  For most bases the user simply supplies the `knots` to be used, which
  must match up with the `k` value supplied (note that the number of
  knots is not always just `k`). See `mgcv` documentation for more
  details.

## Value

a list containing the design matrix `Z`, a (potentially nested) list of
penalty matrices `S`, the `formula`, the `data`, the `knots`, and the
original `mod` object returned by `mgcv`. Note that for tensorproduct
smooths, the corresponding list entry is itself a list, containing the d
marginal penalty matrices if d is the dimension of the tensor product.

## Examples

``` r
data = data.frame(x = runif(100), 
                  y = runif(100),
                  g = factor(rep(1:10, each = 10)))

# unvariate thin plate regression spline
modmat = make_matrices(~ s(x), data)
# univariate P-spline
modmat = make_matrices(~ s(x, bs = "ps"), data)
# adding random intercept
modmat = make_matrices(~ s(g, bs = "re") + s(x, bs = "ps"), data)
# tensorproduct of x and y
modmat = make_matrices(~ s(x) + s(y) + ti(x,y), data)
```
