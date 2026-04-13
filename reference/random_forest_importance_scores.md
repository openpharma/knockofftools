# Internal function called to return the importance scores from random forest

Internal function called to return the importance scores from random
forest

## Usage

``` r
random_forest_importance_scores(X, y, trt, type = "regression", ...)
```

## Arguments

- X:

  original data.frame with "numeric" and "factor" columns only.

- y:

  response vector with `length(y) = nrow(X)`. Accepts "numeric"
  (family="gaussian") or binary "factor" (family="binomial"). Can also
  be a survival object of class Surv as obtained from y =
  survival::Surv(time, status).

- type:

  should be "regression" if y is numeric, "classification" if y is a
  binary factor variable or "survival" if y is a survival object.

- ...:

## Value

importance scores
