# Compute Yule Q coefficient for all gene pairs

Compute Yule Q coefficient for all gene pairs

## Usage

``` r
binary.yule(overlap, mat)
```

## Arguments

- overlap:

  The pairwise overlap matrix

- mat:

  The binary GAM (features x samples)

## Value

Matrix of Yule Q coefficients

## Examples

``` r
am   <- matrix(c(0,1,1,0,1,1), nrow = 2,
               dimnames = list(c("geneA","geneB"), c("s1","s2","s3")))
ovlp <- am.pairwise.alteration.overlap(am)
binary.yule(ovlp, am)
#>       geneA geneB
#> geneA   NaN    -1
#> geneB    -1   NaN
```
