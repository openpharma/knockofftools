# Controls the false discovery rate (FDR) given knockoff W-statistics.

This method was introduced by Candes et al. (2018)

## Usage

``` r
selections_control_FDR(W, level)
```

## Arguments

- W:

  a vector of knockoff W-statistics (feature statistics).

- level:

  The nominal level at which to control the FDR.

## Value

the selected variables

## Details

E. Candes, Y. Fan, L. Janson, & J. Lv, (2018). Panning for
gold:‘model‐X’knockoffs for high dimensional controlled variable
selection. Journal of the Royal Statistical Society: Series B
(Statistical Methodology), 80(3), 551-577.

## Examples

``` r
W <- rnorm(100)
selections_control_FDR(W, level=0.5)
#> integer(0)
```
