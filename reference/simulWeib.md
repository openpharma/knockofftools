# Function that simulates response from Cox model with Weibull baseline hazard:

Function that simulates response from Cox model with Weibull baseline
hazard:

## Usage

``` r
simulWeib(N, lambda0, rho, lp)
```

## Arguments

- N:

  sample size

- lambda0:

  baseline hazard scale parameter

- rho:

  baseline hazard shape parameter

- lp:

  linear predictor

## Value

survival object with simulated event times under mild censoring

## Examples

``` r
# Simulate 10 Gaussian covariate predictors:
X <- generate_X(n=100, p=10, p_b=0, cov_type="cov_equi", rho=0.2)

# create linear predictor with first 5 beta-coefficients = 1 (all other zero)
lp <- generate_lp(X, p_nn = 5, a=1)

# Simulate from Weibull hazard with with baseline hazard h0(t) = lambda*rho*t^(rho-1)
# and linear predictor, whose first 3 coefficients are non-zero:
y <- simulWeib(N=nrow(X), lambda0=0.01, rho=1, lp=lp)
```
