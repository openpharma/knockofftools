# Heuristic check for whether numeric variables can be reasonably treated as continuous

Heuristic check for whether numeric variables can be reasonably treated
as continuous

## Usage

``` r
check_if_continuous(X)
```

## Arguments

- X:

  the design matrix of interest with columns either "numeric" or
  "factor"

## Value

a logical TRUE or FALSE depending on whether n_distinct(x) \> 30
