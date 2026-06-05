# Compute TMB-weighted pairwise alteration overlap

Compute TMB-weighted pairwise alteration overlap

## Usage

``` r
am.weight.pairwise.alteration.overlap(am, W)
```

## Arguments

- am:

  Binary alteration matrix (features x samples)

- W:

  Weight matrix (features x samples) with per-sample TMB weights

## Value

Weighted pairwise overlap matrix (features x features).

## Examples

``` r
am <- matrix(c(0,1,1,0,1,1), nrow = 2,
             dimnames = list(c("geneA","geneB"), c("s1","s2","s3")))
W  <- matrix(1, nrow = 2, ncol = 3,
             dimnames = list(c("geneA","geneB"), c("s1","s2","s3")))
am.weight.pairwise.alteration.overlap(am, W)
#>       geneA geneB
#> geneA     2     1
#> geneB     1     2
```
