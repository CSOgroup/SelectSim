# Initialize an Alteration Matrix Stats (AMS) container

Initialize an Alteration Matrix Stats (AMS) container

## Usage

``` r
new.AMS(am)
```

## Arguments

- am:

  The alteration matrix (checked for NULL)

## Value

Empty AMS list object.

## Examples

``` r
new.AMS(matrix(0, 2, 2))
#> list()
#> attr(,"class")
#> [1] "AMS"
```
