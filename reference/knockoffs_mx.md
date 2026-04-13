# Gaussian MX-knockoffs for continuous variables

Gaussian MX-knockoffs for continuous variables

## Usage

``` r
knockoffs_mx(X)
```

## Arguments

- X:

  data.frame (or tibble) with "numeric" columns only. The number of
  columns, ncol(X) needs to be \> 2.

## Value

Second-order multivariate Gaussian knockoff copy of X

## Details

`knockoffs_mx` performs MX knockoff simulation.

## Examples

``` r
#' library(knockofftools)

set.seed(1)

X <- generate_X(n=100, p=6, p_b=0, cov_type="cov_equi", rho=0.5)

Xk <- knockoffs_mx(X)
```
