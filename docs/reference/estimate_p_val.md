# Compute empirical two-sided p-value for a gene pair

Compute empirical two-sided p-value for a gene pair

## Usage

``` r
estimate_p_val(robs_co, obs.co, gene1, gene2)
```

## Arguments

- robs_co:

  List of null model overlap matrices (one per permutation)

- obs.co:

  Observed pairwise overlap matrix

- gene1:

  Name of the first gene/alteration

- gene2:

  Name of the second gene/alteration

## Value

Two-sided empirical p-value
