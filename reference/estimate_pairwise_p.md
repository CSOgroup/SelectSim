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

## Examples

``` r
# \donttest{
data(luad_run_data, package = "SelectSim")
result <- selectX(M = luad_run_data$M,
                  sample.class = luad_run_data$sample.class,
                  alteration.class = luad_run_data$alteration.class,
                  n.cores = 1, min.freq = 10, n.permut = 10,
                  verbose = FALSE, estimate_pairwise = TRUE)
head(result$result[, c("SFE_1","SFE_2","pairwise_p")])
#>              SFE_1 SFE_2 pairwise_p
#> KRAS - TP53   KRAS  TP53  0.0000000
#> EGFR - KRAS   EGFR  KRAS  0.0000000
#> STK11 - TP53 STK11  TP53  0.0000000
#> KRAS - STK11  KRAS STK11  0.8571429
#> BRAF - KRAS   BRAF  KRAS  0.0000000
#> EGFR - STK11  EGFR STK11  0.0000000
# }
```
