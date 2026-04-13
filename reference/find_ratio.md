# Ratio from Ren et al. (2021) used to derandomized knockoffs.

Do not call this function on its own.

## Usage

``` r
find_ratio(M, eta)
```

## Arguments

- M:

  number of knockoff runs

- eta:

  the proportion of runs in which a feature must be selected to be
  selected in the overall derandomized procedure.

## Value

returns pre-calculated ratios to be used in the multiple knockoff
variable selection with error.type = "pfer" or "kfwer"
