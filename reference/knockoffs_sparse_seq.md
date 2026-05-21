# Sparse sequential knockoff generation algorithm

This function takes as input a data frame X and returns its sparse
sequential knockoff copy. Sparse sequential knockoffs first calculates
the adjacency matrix of X (i.e. identifies the zeros/non-zeros of the
precision matrix of X). Then it proceeds with the usual sequential
knockoffs algorithm, except now each sequential regression only includes
covariates that correspond to non-zero elements of the precision matrix
of X. This reduces the number of covariates per regression in the
original sequential knockoff algorithm. To gain additional speed-up (as
compared to sequential knockoffs) we apply least squares as the default
method for estimation in each regression, unless the number of
covariates exceeds half the number of observations, i.e. p \> n/2. In
this case we apply the elastic net regularized regression.

## Usage

``` r
knockoffs_sparse_seq(X, adjacency.matrix = NULL)
```

## Arguments

- X:

  data.frame (or tibble) with "numeric" and "factor" columns only. The
  number of columns, ncol(X) needs to be \> 2.

- adjacency.matrix:

  optional user specified adjacency matrix (i.e. binary indicator matrix
  corresponding to the non-zero elements of the precision matrix of X).
  Defaults to NULL and is then estimated within the function call.

## Value

sparse sequential knockoff copy of X. A data.frame or tibble of same
type and dimensions as X.

## Examples

``` r
library(knockofftools)

set.seed(1)

X <- generate_X(n=100, p=6, p_b=2, cov_type="cov_equi", rho=0.5)

# knockoffs based on sequential elastic-net regression:
Xk <- knockoffs_sparse_seq(X)
#> # weights:  8 (7 variable)
#> initial  value 69.314718 
#> iter  10 value 42.911908
#> final  value 42.498721 
#> converged
#> # weights:  8 (7 variable)
#> initial  value 69.314718 
#> iter  10 value 48.791354
#> final  value 48.785677 
#> converged
```
