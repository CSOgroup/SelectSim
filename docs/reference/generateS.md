# Generate S matrix

Generate S matrix

## Usage

``` r
generateS(gam, sample.weights, upperBound = 1)
```

## Arguments

- gam:

  the gam with genes\*samples

- sample.weights:

  the samples weights

- upperBound:

  clip the values greater than 1 to keep it bounded between 0 to 1

## Value

S the S matrix

## Examples

``` r
gam <- matrix(c(0,1,1,0,1,1), nrow = 2,
              dimnames = list(c("geneA","geneB"), c("s1","s2","s3")))
generateS(gam, sample.weights = c(s1 = 1, s2 = 1, s3 = 1))
#>              s1        s2        s3
#> geneA 0.6666667 0.6666667 0.6666667
#> geneB 0.6666667 0.6666667 0.6666667
```
