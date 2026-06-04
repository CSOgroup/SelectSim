# Build the full interaction results table from selectX outputs

Build the full interaction results table from selectX outputs

## Usage

``` r
interaction.table(
  al,
  als,
  obs,
  wobs,
  r.obs = NULL,
  r.wobs = NULL,
  null,
  maxFDR = 0.2,
  n.cores = 1,
  estimate_pairwise = FALSE,
  n.permut = 1000
)
```

## Arguments

- al:

  Alteration landscape object (al\$am, etc.)

- als:

  Alteration landscape stats (from al.stats)

- obs:

  Observed pairwise overlap matrix

- wobs:

  Weighted observed pairwise overlap matrix

- r.obs:

  List of null model overlap matrices

- r.wobs:

  List of null model weighted overlap matrices

- null:

  Null model list (simulated binary matrices)

- maxFDR:

  FDR cutoff for calling significant results

- n.cores:

  Number of cores

- estimate_pairwise:

  Whether to compute per-pair empirical p-values

- n.permut:

  Number of permutations

## Value

Data frame with one row per gene pair and columns for overlap, effect
sizes, FDR, and interaction type.
