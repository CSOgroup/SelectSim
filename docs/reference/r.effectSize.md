# Compute effect sizes for null model permutations

Compute effect sizes for null model permutations

## Usage

``` r
r.effectSize(null_overlap, mean_mat, n.permut = 1000, n.cores = 1)
```

## Arguments

- null_overlap:

  List of null model overlap matrices (one per permutation)

- mean_mat:

  Mean overlap matrix across all permutations

- n.permut:

  Number of permutations

- n.cores:

  Number of cores (currently unused; sequential only)

## Value

List of effect size vectors, one per permutation.
