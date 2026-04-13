# Simulate Gaussian and binary covariate predictors

simulate Gaussian predictors with mean zero and covariance structure
determined by "cov_type" argument. Then p_b randomly selected columns
are dichotomized.

## Usage

``` r
generate_X(n, p, p_b, cov_type, rho = 0.5)
```

## Arguments

- n:

  number of observations (rows of X)

- p:

  total number of covariates (columns of X) both continuous and binary

- p_b:

  number of binary covariates (0 \<= p_b \<= p)

- cov_type:

  character string specifying the covariance function. Can be one of
  "cov_diag" (independent columns), "cov_equi" (equi-correlated
  columns), or "cov_ar1" (ar1-correlated columns). The columns are
  shuffled during simulation

- rho:

  correlation parameter; input to the cov_type function

## Value

the simulated data.frame with n rows and p columns (p_b of which are
binary and p-p_b of which are gaussian). Each column is either of class
"numeric" or "factor".

## Details

This function simulates a data frame, whose rows are multivariate
Gaussian with mean zero and covariance structure determined by
"cov_type" argument. Then p_b randomly selected columns are dichotomized
with the function 1(x\>0). The continuous columns are of class "numeric"
and the binary columns are set to class "factor".

## Examples

``` r
library(knockofftools)

# all columns are continuous:
X <- generate_X(n=100, p=6, p_b=0, cov_type="cov_equi", rho=0.5)

round(cor(X), 2)
#>      X1   X2   X3   X4   X5   X6
#> X1 1.00 0.50 0.53 0.59 0.44 0.52
#> X2 0.50 1.00 0.64 0.61 0.53 0.63
#> X3 0.53 0.64 1.00 0.53 0.58 0.61
#> X4 0.59 0.61 0.53 1.00 0.54 0.62
#> X5 0.44 0.53 0.58 0.54 1.00 0.52
#> X6 0.52 0.63 0.61 0.62 0.52 1.00

# two of the six columns are dichotomized (and set to class factor):
X <- generate_X(n=100, p=6, p_b=2, cov_type="cov_equi", rho=0.5)

# The class of each column:
unlist(lapply(X, class))
#>        X1        X2        X3        X4        X5        X6 
#>  "factor" "numeric" "numeric"  "factor" "numeric" "numeric" 
```
