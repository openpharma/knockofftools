# Controls the k-familywise error rate (k-FWER) given a vector of knockoff W-statistics.

This method was introduced by Janson and Su (2016),, while we used the
implementation
https://github.com/zhimeir/derandomized_knockoffs_paper/blob/master/R/vanilla_fwer_filter.R

## Usage

``` r
selections_control_kFWER(W, level, k)
```

## Arguments

- W:

  a vector of knockoff W-statistics (feature statistics)

- level:

  The nominal level at which to control the k-FWER.

- k:

  a positive integer corresponding to k-FWER (multiple testing when one
  seeks to control at least k false discoveries)

## Value

the selected variables

## Details

L. Janson, & W. Su, (2016). Familywise error rate control via knockoffs.
Electronic Journal of Statistics, 10(1), 960-975.

## Examples

``` r
W <- rnorm(100)
selections_control_kFWER(W, level=0.1, k=5)
#>  [1]  2 10 59 30  8 77 58 86 74 54 72 39 49
```
