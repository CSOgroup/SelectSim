# Identify outlier null-model matrices

Flags simulations whose mean per-gene and per-sample deviation from the
observed counts falls in the top 10%, indicating numerical instability.
These are removed before computing effect sizes to prevent them from
inflating the null distribution.

## Usage

``` r
retrieveOutliers(obj, nSim = 1000)
```

## Arguments

- obj:

  SelectX object (list with al and null fields)

- nSim:

  Number of simulations in the null model

## Value

Logical vector; TRUE for each null matrix flagged as an outlier.

## Examples

``` r
# \donttest{
data(luad_run_data, package = "SelectSim")
result <- selectX(M = luad_run_data$M,
                  sample.class = luad_run_data$sample.class,
                  alteration.class = luad_run_data$alteration.class,
                  n.cores = 1, min.freq = 10, n.permut = 10,
                  verbose = FALSE)
outliers <- retrieveOutliers(result$obj, nSim = result$obj$nSim)
sum(outliers)
#> [1] 2
# }
```
