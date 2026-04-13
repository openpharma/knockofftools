# Internal function called to return the importance scores from causal forest

Internal function called to return the importance scores from causal
forest

## Usage

``` r
causal_forest_importance_scores(
  X,
  y,
  trt,
  type = "regression",
  shuffle = FALSE,
  ...
)
```

## Arguments

- X:

  original data.frame with "numeric" and "factor" columns only.

- y:

  response vector with `length(y) = nrow(X)`. Accepts "numeric"
  (family="gaussian") or binary "factor" (family="binomial"). Can also
  be a survival object of class Surv as obtained from y =
  survival::Surv(time, status).

- trt:

  binary treatment indicator

- type:

  should be "regression" if y is numeric, "classification" if y is a
  binary factor variable or "survival" if y is a survival object.

- shuffle:

  boolean variable, that takes the value FALSE if we do not want the
  target y to be shuffled, and TRUE otherwise

- ...:

## Value

importance scores
