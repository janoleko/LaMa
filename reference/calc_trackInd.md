# Calculate the index of the first observation of each track based on an ID variable

Function to conveniently calculate the trackInd variable that is needed
internally when fitting a model to longitudinal data with multiple
tracks.

## Usage

``` r
calc_trackInd(ID)
```

## Arguments

- ID:

  ID variable of track IDs that is of the same length as the data to be
  analysed

## Value

A vector of indices of the first observation of each track which can be
passed to the forward and forward_g to sum likelihood contributions of
each track

## Examples

``` r
uniqueID = c("Animal1", "Animal2", "Animal3")
ID = rep(uniqueID, c(100, 200, 300))
trackInd = calc_trackInd(ID)
```
