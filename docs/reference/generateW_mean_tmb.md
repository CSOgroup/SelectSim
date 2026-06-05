# Generate sample weight matrix from TMB values

Computes a per-sample weight matrix based on the ratio of each sample's
TMB to the expected (mean) TMB. Samples with higher-than-expected TMB
receive lower weights via a penalty `lambda`, controlled by the
fold-change threshold `tau`.

## Usage

``` r
generateW_mean_tmb(
  tmb,
  mean_tmb,
  ngenes,
  lambda = 0.3,
  tau = 1,
  discrete = TRUE
)
```

## Arguments

- tmb:

  Numeric vector of per-sample TMB values.

- mean_tmb:

  Numeric scalar; the reference (expected) TMB used to compute fold
  changes.

- ngenes:

  Integer; number of genes (rows) in the output weight matrix.

- lambda:

  Numeric; weight penalty factor. Higher values penalise high-TMB
  samples more strongly (default 0.3).

- tau:

  Numeric; fold-change threshold below which no penalty is applied
  (default 1).

- discrete:

  Logical; if `TRUE`, fold changes are rounded up before applying the
  penalty (default `TRUE`).

## Value

Numeric matrix of sample weights (ngenes x length(tmb)).

## Examples

``` r
tmb      <- c(s1 = 10, s2 = 50, s3 = 20)
mean_tmb <- 25
generateW_mean_tmb(tmb, mean_tmb, ngenes = 3)
#>      [,1]      [,2] [,3]
#> [1,]    1 0.7692308    1
#> [2,]    1 0.7692308    1
#> [3,]    1 0.7692308    1
```
