# Identify outlier null-model matrices

Flags simulations whose mean per-gene and per-sample deviation from the
observed counts falls in the top 10%, indicating numerical instability.

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
