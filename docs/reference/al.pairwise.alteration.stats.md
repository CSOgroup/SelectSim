# Compute pairwise alteration statistics for an alteration landscape

Compute pairwise alteration statistics for an alteration landscape

## Usage

``` r
al.pairwise.alteration.stats(al, als = NULL, do.blocks = FALSE)
```

## Arguments

- al:

  SelectX object (list containing al, W, etc.)

- als:

  Alteration landscape stats (from al.stats); computed internally if
  NULL

- do.blocks:

  Whether to also compute block-level pairwise stats

## Value

List with overlap and w_overlap matrices, plus optional sample.blocks
entries.
