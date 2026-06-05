# Compute summary statistics for a binary alteration matrix

Compute summary statistics for a binary alteration matrix

## Usage

``` r
am.stats(am)
```

## Arguments

- am:

  Binary alteration matrix (features x samples)

## Value

AMS object with basic counts: n.samples, n.alterations, n.occurrences,
per-sample and per-feature counts.

## Examples

``` r
am <- matrix(c(0,1,1,0,1,1), nrow = 2,
             dimnames = list(c("geneA","geneB"), c("s1","s2","s3")))
am.stats(am)
#> $n.samples
#> [1] 3
#> 
#> $n.alterations
#> [1] 2
#> 
#> $n.occurrences
#> [1] 4
#> 
#> $alterations.per.sample
#> s1 s2 s3 
#>  1  1  2 
#> 
#> $alteration.count
#> geneA geneB 
#>     2     2 
#> 
#> attr(,"class")
#> [1] "AMS"
```
