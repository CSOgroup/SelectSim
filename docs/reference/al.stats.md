# Compute alteration landscape statistics

Compute alteration landscape statistics

## Usage

``` r
al.stats(al)
```

## Arguments

- al:

  SelectX object (list containing al, W, etc. as returned by selectX)

## Value

ALS object with overall and per-block alteration statistics.

## Examples

``` r
# \donttest{
data(luad_run_data, package = "SelectSim")
result <- selectX(M = luad_run_data$M,
                  sample.class = luad_run_data$sample.class,
                  alteration.class = luad_run_data$alteration.class,
                  n.cores = 1, min.freq = 10, n.permut = 10,
                  verbose = FALSE)
stats <- al.stats(result$obj)
# }
```
