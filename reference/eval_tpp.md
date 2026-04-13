# True positive proportion (tpp) as function of selection and known positives:

True positive proportion (tpp) as function of selection and known
positives:

## Usage

``` r
eval_tpp(selected, positives)
```

## Arguments

- selected:

  vector of indices of selected variables

- negatives:

  vector of indices of known non-null variables (that influence
  response)

## Value

true positive rate

## Examples

``` r
#' library(knockofftools)

eval_tpp(selected=c(1,2,3,5), positives=1:4)
#> [1] 0.75
```
