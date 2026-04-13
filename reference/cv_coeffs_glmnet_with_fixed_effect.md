# Internal function called within the stat_glmnet.

Do not call this function on its own. Fits cross-validated glmnet model
with fixed effect.

## Usage

``` r
cv_coeffs_glmnet_with_fixed_effect(
  X_fixed,
  X,
  y,
  family,
  nlambda = 500,
  penalty.factor,
  ...
)
```

## Arguments

- X:

  original data.frame (or tibble) with "numeric" and "factor" columns
  only. The number of columns, ncol(X) needs to be \> 2.

- y:

  response vector with `length(y) = nrow(X)`. Accepts "numeric"
  (family="gaussian") or binary "factor" (family="binomial"). Can also
  be a survival object of class Surv as obtained from y =
  survival::Surv(time, status).

- family:

  should be "gaussian" if y is numeric, "binomial" if y is a binary
  factor variable or "cox" if y is a survival object.

- nlambda:

  length of lambda penalty sequence

- penalty.factor:

  Separate penalty factors, passed to the glmnet::cv.glmnet function,
  can be applied to each coefficient. This is a number that multiplies
  lambda to allow differential shrinkage. Can be 0 for some variables,
  which implies no shrinkage, and that variable is always included in
  the model.

- ...:

  other parameters passed to glmnet::cv.glmnet

- X.fixed:

  a data.frame (or tibble) with "numeric" and "factor" columns
  corresponding to covariates or terms that should be treated as fixed
  effects in the model.

## Value

coefficients of both the original and knockoff variables
