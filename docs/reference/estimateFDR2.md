# Estimate FDR by scanning observed vs null effect sizes

Estimate FDR by scanning observed vs null effect sizes

## Usage

``` r
estimateFDR2(obs, exp, nSim, maxFDR = 0.25)
```

## Arguments

- obs:

  Vector of observed effect sizes

- exp:

  Vector of null model effect sizes (all permutations concatenated)

- nSim:

  Number of permutations used to generate exp

- maxFDR:

  FDR cutoff; scanning stops once FDR exceeds this value

## Value

Vector of FDR values, one per element of obs.

## Examples

``` r
set.seed(1)
obs <- c(0.8, 0.5, 0.3, 0.1)
exp <- runif(400)
estimateFDR2(obs, exp, nSim = 100)
#> [1] 0.72 1.00 1.00 1.00
```
