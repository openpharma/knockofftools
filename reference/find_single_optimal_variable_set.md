# Select variables based on (heuristic) mode of multiple variable selections

Do not call this function on its own

## Usage

``` r
find_single_optimal_variable_set(S, p, trim = 0.5)
```

## Arguments

- S:

  list of variable selection indices

- p:

  number of variables. Each element of the list of selection indices
  should be a subset of 1:p.

- trim:

  trimming probability threshold. A sensible default is `trim=0.5`.

## Value

a single "most frequent" variable selection among the multiple
selections in S.
