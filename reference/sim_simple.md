# Simple knockoff generator based on least squares fit (continuous variables) or multinomial logistic regression (factor variables) respectively. If X is empty, knockoffs are sampled from the marginal distribution of y

Simple knockoff generator based on least squares fit (continuous
variables) or multinomial logistic regression (factor variables)
respectively. If X is empty, knockoffs are sampled from the marginal
distribution of y

## Usage

``` r
sim_simple(y, X)
```

## Arguments

- y:

  response vector (either "numeric" or "factor") that gets passed to
  cv.glmnet

- X:

  data.frame of covariates that are passed to cv.glmnet

## Value

simulated response
