# Scatter plot of observed vs expected weighted co-mutation

Plots each gene pair as a point with observed weighted co-mutation on
the y-axis and expected (null model mean) on the x-axis. Significant
co-mutations (CO) and mutual exclusivities (ME) are coloured;
non-significant pairs are grey.

## Usage

``` r
obs_exp_scatter(result, title)
```

## Arguments

- result:

  Result table from `selectX()$result`.

- title:

  Plot title.

## Value

A ggplot2 object.

## Examples

``` r
# \donttest{
data(luad_result, package = "SelectSim")
obs_exp_scatter(result = luad_result, title = "TCGA LUAD")

# }
```
