# Knockoff (feature) statistics for a single knockoff

Do not call this function on its own. Use `knockoff.statistics` instead.

## Usage

``` r
.knockoff.statistics.single(
  y,
  X,
  type = "regression",
  knockoff.method = "seq",
  statistic = "stat_glmnet",
  trt = NULL,
  adjacency.matrix = NULL,
  ...
)
```

## Arguments

- y:

  response vector with `length(y) = nrow(X)`. Accepts "numeric", binary
  "factor", or survival ("Surv") object.

- X:

  data.frame (or tibble) with "numeric" and "factor" columns only. The
  number of columns, ncol(X) needs to be \> 2.

- type:

  should be "regression" if y is numeric, "classification" if y is a
  binary factor variable or "survival" if y is a survival object.

- knockoff.method:

  what type of knockoffs to calculate. Defaults to sequential knockoffs,
  knockoff.method="seq", but other options are "sparseseq" and "mx". The
  "mx" option only works if all columns of X are continuous.

- statistic:

  knockoff feature statistic function, defaults to glmnet coefficient
  difference (statistic="stat_glmnet"; see ?stat_glmnet). Other options
  include statistic="stat_random_forest" (see ?stat_random_forest),
  statistic="stat_predictive_glmnet" (see ?stat_predictive_glmnet) or
  statistic="stat_predictive_causal_forest" (see
  ?stat_predictive_causal_forest).

- trt:

  binary treatment (factor) variable required if statistic involves a
  predictive knockoff filter (i.e. if statistic="stat_predictive_glmnet"
  or statistic="stat_predictive_causal_forest")

- adjacency.matrix:

  optional user specified adjacency matrix (i.e. binary indicator matrix
  corresponding to the non-zero elements of the precision matrix of X).
  Defaults to NULL and is then estimated within the function call.

- ...:

  additional parameters passed to the "statistic" function (note that
  the knockoffs parameter X_k should not be entered by user; it is
  already calculated inside the knockoff.statistics function).

## Value

data.frame with a single knockoff statistics W as column.
