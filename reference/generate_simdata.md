# Generate simulated data set for vignettes

Generate simulated data set for vignettes

## Usage

``` r
generate_simdata()
```

## Value

generates the simulated data as obtained by the call data(simdata)

## Examples

``` r
library(knockofftools)
data(simdata)
all.equal(generate_simdata(), simdata)
#> [1] TRUE
```
