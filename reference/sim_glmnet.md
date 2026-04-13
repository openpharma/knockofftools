# Simulate from glmnet penalized regression model

Simulate from glmnet penalized regression model

## Usage

``` r
sim_glmnet(y, X, ...)
```

## Arguments

- y:

  response vector (either "numeric" or "factor") that gets passed to
  cv.glmnet

- X:

  data.frame of covariates that are passed to cv.glmnet

- ...:

  other parameters passed to the function cv.glmnet

## Value

simulated response

## Examples

``` r
library(knockofftools)

set.seed(1)

X = data.frame(matrix(rnorm(100 * 20), 100, 20))
y = X[,1] + rnorm(100)

# simulate from elastic-net regression:
ysim = sim_glmnet(y=y, X=X)

# simulated versus input response:
plot(y, ysim)
```
