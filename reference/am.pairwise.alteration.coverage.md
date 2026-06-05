# Compute pairwise alteration coverage statistics

Compute pairwise alteration coverage statistics

## Usage

``` r
am.pairwise.alteration.coverage(overlap_M, M.stats, w_overlap_M)
```

## Arguments

- overlap_M:

  The pairwise overlap matrix

- M.stats:

  The alteration matrix stats (from am.stats)

- w_overlap_M:

  The weighted pairwise overlap matrix

## Value

List with overlap and w_overlap sparse matrices.

## Examples

``` r
am <- matrix(c(0,1,1,0,1,1), nrow = 2,
             dimnames = list(c("geneA","geneB"), c("s1","s2","s3")))
W  <- matrix(1, nrow = 2, ncol = 3,
             dimnames = list(c("geneA","geneB"), c("s1","s2","s3")))
ovlp  <- am.pairwise.alteration.overlap(am)
wovlp <- am.weight.pairwise.alteration.overlap(am, W)
stats <- am.stats(am)
am.pairwise.alteration.coverage(ovlp, stats, wovlp)
#> $overlap
#> 2 x 2 Matrix of class "dsyMatrix"
#>       geneA geneB
#> geneA     2     1
#> geneB     1     2
#> 
#> $w_overlap
#> 2 x 2 Matrix of class "dsyMatrix"
#>       geneA geneB
#> geneA     2     1
#> geneB     1     2
#> 
```
