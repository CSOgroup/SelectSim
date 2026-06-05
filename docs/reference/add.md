# Sum a list of matrices element-wise

Sum a list of matrices element-wise

## Usage

``` r
add(x)
```

## Arguments

- x:

  List of matrices of identical dimensions

## Value

Single matrix that is the element-wise sum of all matrices in x

## Examples

``` r
mats <- list(matrix(1:4, 2, 2), matrix(1:4, 2, 2))
add(mats)
#>      [,1] [,2]
#> [1,]    2    6
#> [2,]    4    8
```
