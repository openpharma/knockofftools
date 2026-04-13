# Knockoffs for general covariate data frames

Knockoffs for general covariate data frames

## Usage

``` r
knockoff(X, method = "seq", ...)
```

## Arguments

- X:

  data.frame (or tibble) with "numeric" and "factor" columns only. The
  number of columns, ncol(X) needs to be \> 2.

- method:

  what type of knockoffs to calculate. Defaults to sequential knockoffs,
  method="seq", which works for both numeric and factor variables. The
  other options are method="sparseseq", which is the sparse sequential
  knockoff generation algorithm, and method="mx", which only works if
  all columns of X are continuous.

- ...:

  additional parameters passed to the downstream knockoffs generating
  functions

## Value

knockoff copy of X

## Examples

``` r
library(knockofftools)

set.seed(1)

X <- generate_X(n=100, p=6, p_b=2, cov_type="cov_equi", rho=0.5)

# sequential knockoff:
Xk <- knockoff(X)

X <- generate_X(n=100, p=6, p_b=0, cov_type="cov_equi", rho=0.5)

# MX-knockoff:
Xk <- knockoff(X, method="mx")
```
