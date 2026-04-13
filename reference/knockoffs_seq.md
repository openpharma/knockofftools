# Sequential knockoffs for continuous and categorical variables

Sequential knockoffs for continuous and categorical variables

## Usage

``` r
knockoffs_seq(X, seq_simulator = sim_glmnet, ...)
```

## Arguments

- X:

  data.frame (or tibble) with "numeric" and "factor" columns only. The
  number of columns, ncol(X) needs to be \> 2.

- seq_simulator:

  function that simulates sequential knockoffs. Default is the function
  `sim_glmnet`, which simulates response from an estimated elastic-net
  model

- ...:

  other parameters passed to the function seq_simulator. For the default
  (elastic-net sequential seq_simulator, `seq_simulator = sim_glmnet`)
  these other parameters are passed to cv.glmnet.

## Value

sequential knockoff copy of X. A data.frame or tibble of same type and
dimensions as X.

## Details

`knockoffs_seq` performs sequential knockoff simulation using
elastic-net regression.

## Examples

``` r
library(knockofftools)

set.seed(1)

X <- generate_X(n=100, p=6, p_b=2, cov_type="cov_equi", rho=0.5)

# knockoffs based on sequential elastic-net regression:
Xk <- knockoffs_seq(X)
```
