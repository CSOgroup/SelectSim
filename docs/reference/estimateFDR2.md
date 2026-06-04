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
