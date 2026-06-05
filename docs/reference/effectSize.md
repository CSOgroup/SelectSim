# Compute effect size between observed and expected overlap

Compute effect size between observed and expected overlap

## Usage

``` r
effectSize(obs, exp)
```

## Arguments

- obs:

  The observed overlap values

- exp:

  The expected (null model mean) overlap values

## Value

Effect size value(s)

## Examples

``` r
effectSize(obs = 5, exp = 3)
#> [1] 1.414214
effectSize(obs = c(5, 2, 0), exp = c(3, 3, 1))
#> [1]  1.4142136 -0.7071068 -0.7071068
```
