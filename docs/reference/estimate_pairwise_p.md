# Compute p-values for all gene pairs in a results table

Compute p-values for all gene pairs in a results table

## Usage

``` r
estimate_pairwise_p(obs, exp, results, nSim)
```

## Arguments

- obs:

  Observed pairwise overlap matrix

- exp:

  List of null model overlap matrices (one per permutation)

- results:

  Results data frame with SFE_1 and SFE_2 columns

- nSim:

  Number of permutations

## Value

Vector of p-values, one per row in results.
