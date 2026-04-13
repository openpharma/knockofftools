# Simulate Gaussian response from a sparse regression model

Simulate Gaussian response from a sparse regression model

## Usage

``` r
generate_y(X, p_nn, a)
```

## Arguments

- X:

  data.frame with numeric and factor columns only.

- p_nn:

  number of non-null covariate predictors. The regression coefficients
  (beta) corresponding to columns 1:p_nn of x will be non-zero, all
  other are set to zero.

- a:

  amplitude of non-null regression coefficients

## Value

simulated Gaussian response from regression model y = x x is the
(scaled) model.matrix of X.

## Details

This function takes as input data.frame X (created with the function
`generate_X`) that may consist of both numeric and binary factor
columns. This data frame is then expanded to a model matrix x (with the
model.matrix function) and subsequently scaled in the same way as LASSO
scaling. Next we simulate y ~ N(x the remaining coefficients
(p_nn+1):ncol(x) are set to zero.

## Examples

``` r
library(knockofftools)

set.seed(1)

# Simulate 4 Gaussian and 2 binary covariate predictors:
X <- generate_X(n=100, p=6, p_b=2, cov_type="cov_equi", rho=0.5)

# Simulate response from model y = 2*X[,1] + 2*X[,2] + epsilon, where epsilon ~ N(0,1)
y <- generate_y(X, p_nn=2, a=2)
```
