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

## Examples

``` r
set.seed(1)
genes <- c("geneA", "geneB")
robs_co <- lapply(seq_len(20), function(i) {
  m <- matrix(sample(0:5, 4, replace = TRUE), 2, 2,
              dimnames = list(genes, genes))
  m
})
obs_co <- matrix(c(5, 3, 3, 4), 2, 2, dimnames = list(genes, genes))
estimate_p_val(robs_co, obs_co, "geneA", "geneB")
#> [1] 0.8
```
