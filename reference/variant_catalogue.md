# OncoKB v3.9 cancer genes

A dataframe cancer genes with missense mutation annotations

## Usage

``` r
data(variant_catalogue)
```

## Format

A data frame

## Value

A data frame containing cancer-gene and variant annotations.

## Examples

``` r
data(variant_catalogue)
dim(variant_catalogue)
#> [1] 2478    3
head(variant_catalogue)
#>   gene  mut  oncogenic
#> 1 ABL1 M244 Resistance
#> 2 ABL1 L248 Resistance
#> 3 ABL1 G250 Resistance
#> 4 ABL1 Q252 Resistance
#> 5 ABL1 Y253 Resistance
#> 7 ABL1 E255 Resistance
```
