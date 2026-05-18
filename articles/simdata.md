# The simdata dataset

## Introduction

The knockofftools package comes with a simulated data set called
`simdata`. This short vignette demonstrates how the data set was
generated.

``` r

library(knockofftools)
```

## Data generation with `generate_simdata`

The function `generate_simdata` was used to generate the data set
`simdata` that comes with the R-package.

``` r

generate_simdata <- function() {

  RNGkind("L'Ecuyer-CMRG")
  set.seed(56969)

  N <- 2000
  p <- 30
  p_b = 10
  p_nn <- 10

  # Generate a 2000 x 30 Gaussian data.frame under equi-correlation(rho=0.5) structure,
  # with 10 of the columns dichotomized
  X <- generate_X(n=N, p=p, p_b=p_b, cov_type = "cov_equi", rho=0.5)

  # Generate linear predictor lp = X%*%beta where first 10 beta-coefficients are = a, all other = 0.
  lp <- generate_lp(X=X, p_nn=p_nn, a=1)

  # simulate gaussian response with mean = lp and sd = 1.
  Yg <- lp + rnorm(N)

  # simulate bernoulli response with mean = exp(lp)/(1+exp(lp)).
  Yb <- factor(rbinom(N, size=1, prob=exp(lp)/(1+exp(lp))))

  # simulate censored survival times from Cox regression with linear predictor lp:
  Tc <- simulWeib(N=N, lambda = 0.01, rho = 1, lp = lp)

  dat <- data.frame(Yg, Yb, Tc, X)

  return(dat)

}
```

This function simulates a toy data set with 30 covariates
$`X_1, \dots, X_{30}`$, one continuous response $`y_g`$, one binary
response $`y_b`$ and one set of censored survival times $`T_c`$. Let’s
now go through the individual components of the function one by one.

*Simulation of $`X`$:* The `generate_X` function simulates the rows of
an $`n \times p`$ data frame $`X`$ independently from a multivariate
Gaussian distribution with mean $`0`$ and $`p \times p`$ covariance
matrix
``` math
\Sigma_{ij} = \left \{
\begin{array}{lr}
1\{i = j\}, & \text{Independent,} \\
\rho^{1\{i \neq j\}}, & \text{Equicorrelated,} \\
\rho^{|i-j|}, & \text{AR1},
\end{array}
\right.
```
where $`p_b`$ randomly selected columns are then dichotomized with the
indicator function $`\delta(x)=1(x > 0)`$.

``` r

X <- generate_X(n=N, p=p, p_b=p_b, cov_type = "cov_equi", rho=0.5)
```

The covariance type is specified with the parameter `cov_type` and the
correlation coefficient with `rho`. Each column of the resulting
data.frame is either of class `"numeric"` (for the continuous columns)
or `"factor"` (for the binary columns).

*Calculation of linear predictor:* The `generate_lp` function calculates
the linear predictor $`\ell_p=X\beta`$ under sparsity, where the first
`p_nn` regression coefficients are non-zero, all other are set to zero.
The (common) amplitude of the non-zero regression coefficients is
specified with `a`. Here we generate $`\ell_p`$ that implies association
with the first 10 covariates, each with amplitude $`a=1`$.

``` r

lp <- generate_lp(X=X, p_nn=p_nn, a=1)
```

Note that inside `generate_lp` the model.matrix of `X` is first scaled.

*Simulation of Gaussian response:* $`Y_g`$ is Gaussian with mean
$`\mu = \ell_p`$ and standard deviation $`\sigma = 1`$:

``` r

Yg <- lp + rnorm(N)
```

*Simulation of Bernoulli response:* $`Y_b`$ is Bernoulli with success
probability $`\mu = exp(\ell_p)/(1+exp(\ell_p))`$:

``` r

Yb <- factor(rbinom(N, size=1, prob=exp(lp)/(1+exp(lp))))
```

*Simulation of event and censoring times:* The final command:

``` r

Tc <- simulWeib(N=N, lambda = 0.01, rho = 1, lp = lp)
```

generates censored survival times $`T_c`$ from a Cox regression model
``` math
\lambda(t) = \lambda_0(t) \exp \left(\ell_p \right),
```
with Weibull baseline hazard:
``` math
\lambda_0(t) = \lambda_0 \rho t^{\rho-1}
```
where $`\lambda_0 > 0`$ and $`\rho > 0`$ are scale and shape parameters,
respectively. The censoring times $`C`$ are randomly drawn from an
exponential distribution with a small (fixed) rate $`\lambda_C=0.0005`$,
which results in very mild censoring. Once $`T`$ and $`C`$ have been
simulated the function returns a survival object (`Surv`) with
`time = min(T, C)` and `event = 1{T < C}`:

``` r

simulWeib <- function(N, lambda0, rho, lp) {
    lambdaC = 5e-04
    v <- runif(n = N)
    Tlat <- (-log(v)/(lambda0 * exp(lp)))^(1/rho)
    C <- rexp(n = N, rate = lambdaC)
    time <- pmin(Tlat, C)
    status <- as.numeric(Tlat <= C)
    survival::Surv(time = time, event = status)
}
```

## The data set `simdata`

Now let’s have a look at the first few columns of the data set `simdata`

``` r

data(simdata)
```

``` r

head(simdata[,1:9])
#>          Yg Yb     Tc.time   Tc.status X1          X2 X3          X4 X5 X6
#> 1 -8.616078  0 2961.163522    0.000000  0 -1.34958449  1 -2.06459293  0  0
#> 2  2.145781  1   14.818694    1.000000  0  1.17359743  1  0.84017766  1  0
#> 3  2.855925  1    5.662159    1.000000  0  0.06617883  1  0.56828582  1  0
#> 4 -7.393736  0  780.165789    0.000000  0  0.60006023  0 -0.88344294  0  0
#> 5 -5.357337  0 1290.810345    0.000000  0  0.68195033  0 -0.07230342  0  1
#> 6  3.720681  1    3.590458    1.000000  1 -0.65766424  1  0.46974365  0  1
```

and finally confirm that
[`generate_simdata()`](https://openpharma.github.io/knockofftools/reference/generate_simdata.md)
indeed reproduces the `simdata` dataset:

``` r

all.equal(simdata, generate_simdata())
#> [1] TRUE
```
