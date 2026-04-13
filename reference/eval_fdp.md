# False discovery proportion (fdp) as function of selection and known negatives:

False discovery proportion (fdp) as function of selection and known
negatives:

## Usage

``` r
eval_fdp(selected, negatives)
```

## Arguments

- selected:

  vector of indices of selected variables

- negatives:

  vector of indices of known null variables (that don't influence
  response)

## Value

false discovery rate

## Examples

``` r
library(knockofftools)

eval_fdp(selected=c(1,2,3,5), negatives=5:10)
#> [1] 0.25
```
