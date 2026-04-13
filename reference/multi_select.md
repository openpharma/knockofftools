# Select variables based on the heuristic multiple selection algorithm from Kormaksson et al. 'Sequential knockoffs for continuous and categorical predictors: With application to a large psoriatic arthritis clinical trial pool.' Statistics in Medicine. 2021;1–16.

Select variables based on the heuristic multiple selection algorithm
from Kormaksson et al. 'Sequential knockoffs for continuous and
categorical predictors: With application to a large psoriatic arthritis
clinical trial pool.' Statistics in Medicine. 2021;1–16.

## Usage

``` r
multi_select(S, trim = 0.5)
```

## Arguments

- S:

  the binary matrix of selections

- trim:

  trimming probability threshold. A sensible default is `trim=0.5`.

## Value

a single "most frequent" variable selection among the multiple
selections in S.

## Details

M. Kormaksson, L. J. Kelly, X. Zhu, S. Haemmerle, L. Pricop, & D.
Ohlssen (2021). Sequential knockoffs for continuous and categorical
predictors: With application to a large psoriatic arthritis clinical
trial pool. Statistics in Medicine, 40(14), 3313-3328.

## Examples

``` r
library(knockofftools)

set.seed(1)

p = 31
Nknockoff = 100
S <- matrix(sample(0:1,p*Nknockoff, replace=TRUE), p, Nknockoff)

multi_select(S)
#> [1]  1  8  9 13 14 15 21 26 28
```
