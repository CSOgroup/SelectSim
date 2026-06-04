# Compute TMB-weighted pairwise alteration overlap

Compute TMB-weighted pairwise alteration overlap

## Usage

``` r
am.weight.pairwise.alteration.overlap(am, W)
```

## Arguments

- am:

  Binary alteration matrix (features x samples)

- W:

  Weight matrix (features x samples) with per-sample TMB weights

## Value

Weighted pairwise overlap matrix (features x features).
