# Compute pairwise alteration co-occurrence counts

Compute pairwise alteration co-occurrence counts

## Usage

``` r
am.pairwise.alteration.overlap(am)
```

## Arguments

- am:

  Binary alteration matrix (features x samples)

## Value

Square matrix of pairwise co-occurrence counts (features x features).

## Examples

``` r
am <- matrix(c(0,1,1,0,1,1), nrow = 2,
             dimnames = list(c("geneA","geneB"), c("s1","s2","s3")))
am.pairwise.alteration.overlap(am)
#>       geneA geneB
#> geneA     2     1
#> geneB     1     2
```
