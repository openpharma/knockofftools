# Heuristic check for whether a variable can be reasonably treated as continuous

Heuristic check for whether a variable can be reasonably treated as
continuous

## Usage

``` r
check_if_continuous(X)
```

## Arguments

- x:

  a numeric variable vector

## Value

a logical TRUE or FALSE depending on whether n_distinct(x) \> 30
