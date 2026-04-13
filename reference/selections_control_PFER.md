# Controls the per-familywise error rate (PFER) given knockoff W-statistics.

This method was introduced by Janson and Su (2016), while we used the
implementation of
https://github.com/zhimeir/derandomized_knockoffs_paper/blob/master/R/pfer_filter.R

## Usage

``` r
selections_control_PFER(W, level)
```

## Arguments

- W:

  a vector of knockoff W-statistics (feature statistics)

- level:

  The nominal level at which to control the PFER.

## Value

the selected variables

## Details

L. Janson, & W. Su, (2016). Familywise error rate control via knockoffs.
Electronic Journal of Statistics, 10(1), 960-975.

## Examples

``` r
W <- rnorm(100)
selections_control_PFER(W, level=2)
#> [1] 91 71 54
```
