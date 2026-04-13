# Deriving predictive variables using knockoffs

## Introduction

In this vignette we demonstrate how we can use knockofftools to identify
variables that modify the treatment effect. In particular firstly we
will introduce what predictive variables are, and then we will present
various knockoff based filter tailored to the task identifying
predictive variables. More details on this topic can be found in
([Sechidis, Kormaksson, and Ohlssen 2021](#ref-predknockoffpaper))

When we perform an intervention there are typically two types of
variables we are interested in. First, prognostic variables, which,
independent of the treatment, have a direct impact on the target $Y.$
Second, predictive variables, which have an impact on the outcome by
interacting with the treatment and explain the heterogeneous treatment
effect. To illustrate the distinction between these two types of
variables, a common approach is to assume a known data generating model
([Lipkovich, Dmitrienko, and B D’Agostino Sr
2017](#ref-lipkovich2017tutorial)): $$\begin{array}{r}
{{\mathbb{E}}\left( Y|X = \mathbf{x},T = t \right) = h(\mathbf{x}) + g(\mathbf{x})t.}
\end{array}$$ Prognostic variables are those that contribute to
$h(\mathbf{x})$ (i.e., ‘main effects’), while predictive are the ones
that contribute to $g(\mathbf{x})$. The focus of this vignette will be
on the predictive variables. Here we will show how to use use two
different knockoff filters to discover predictive variables, while
controlling type-I error.

## Generate synthetic data

Before showing how to use the fitlers we will generate some synthetic
data, where we will know thw ground truth, i.e. which variables are
predictive. Firstly, we will load the package

``` r
library(knockofftools)
```

Then we will simulate 50 gaussian covariate predictors.

``` r
X <- generate_X(n=500, p=10, p_b=0, cov_type="cov_diag", rho=0.2)
```

Then we will calculate the linear predictor for the prognostic part
($h(\mathbf{x})$), where the first `p_nn` regression coefficients are
non-zero, all other are set to zero. The (common) amplitude of the
non-zero regression coefficients is specified with `a`. In other words
the first `p_nn` variables are the *prognostic variables*.

``` r
lp <- generate_lp(X, p_nn = 5, a=1)
```

Then we will generate a binary variable that will serve as the treatment
indicator.

``` r
trt = sample(c(1,0), nrow(X), replace=TRUE)
```

We will add in the linear predictor some interaction terms, in this case
the 6th and 7th variable are the *predictive variables*.

``` r
lp.pred = lp + 1*trt*( as.integer(X[,6]>0) + as.integer(X[,7]>0))
```

And finally, we will generate the outcome variable

``` r
y <- lp.pred + rnorm(nrow(X))
```

So we want our filters to return the 6th and 7th variable, since are the
predictive ones, and at the same time provide the nominal control.

## Filter 1: Using LASSO regression coefficients of the treatment interaction terms

One direct way of deriving importance scores that capture the predictive
strength of each variable is to use a LASSO linear model, modelling both
main and interaction effects, and check the absolute value of the
coefficient of the interaction terms. In `knockofftools` the function
`stat_predictive_glmnet` implements the filter that follows this
approach.

``` r
W <- knockoff.statistics(y=y, X=X, type="regression", statistic = "stat_predictive_glmnet", trt=trt, M=10)
#> Running sequentially ('LOCAL') ...
S = variable.selections(W, level = 1, error.type="pfer")
S$stable.variables
#> [1] "X4" "X6" "X7"
```

## Filter 2: Using importance scores derived from causal forest

Causal trees (CT), introduced by Athey and Imbens ([Athey and Imbens
2016](#ref-athey2016recursive)), is a recursive partitioning algorithm
for estimating heterogeneous treatment effects. The main difference
between CT and Classification and regression trees (CART) is that CART
focus on predicting the outcome $Y,$ while CT focus on estimating the
treatment effect. % for each example (patient). In CT, heterogeneous
effects are described by the conditional average treatment effect
(CATE):
${\mathbb{E}}\left\lbrack Y(1) - Y(0)|X = \mathbf{x} \right\rbrack$.

We have implemented a knockoff filter that uses the importance scores
from the causal forest to derive predictive markers.

``` r
W <- knockoff.statistics(y=y, X=X, type="regression", statistic = "stat_predictive_causal_forest", trt=trt, M=10)
#> Running sequentially ('LOCAL') ...
S = variable.selections(W, level = 1, error.type="pfer")
S$stable.variables
#> [1] "X7"
```

## References

Athey, Susan, and Guido Imbens. 2016. “Recursive Partitioning for
Heterogeneous Causal Effects.” *Proceedings of the National Academy of
Sciences* 113 (27): 7353–60.

Lipkovich, Ilya, Alex Dmitrienko, and Ralph B D’Agostino Sr. 2017.
“Tutorial in Biostatistics: Data-Driven Subgroup Identification and
Analysis in Clinical Trials.” *Statistics in Medicine* 36 (1): 136–96.

Sechidis, K., M. Kormaksson, and D. Ohlssen. 2021. “Using Knockoffs for
Controlled Predictive Biomarker Identification.” *Statistics in
Medicine*.
