# Generate linear predictor with first p_nn beta coefficients = a, all other = 0

Generate linear predictor with first p_nn beta coefficients = a, all
other = 0

## Usage

``` r
generate_lp(X, p_nn, a)
```

## Arguments

- X:

  data.frame with numeric and factor columns only

- p_nn:

  number of non-null covariate predictors. The regression coefficients
  (beta) corresponding to columns 1:p_nn of x will be non-zero, all
  other are set to zero.

- a:

  amplitude of non-null regression coefficients

## Value

linear predictor X
