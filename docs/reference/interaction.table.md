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

## Examples

``` r
# \donttest{
data(luad_run_data, package = "SelectSim")
result <- selectX(M = luad_run_data$M,
                  sample.class = luad_run_data$sample.class,
                  alteration.class = luad_run_data$alteration.class,
                  n.cores = 1, min.freq = 10, n.permut = 10,
                  verbose = FALSE)
head(result$result)
#>              SFE_1 SFE_2         name support_1 support_2     freq_1    freq_2
#> KRAS - TP53   KRAS  TP53  KRAS - TP53       154       221 0.30677291 0.4402390
#> EGFR - KRAS   EGFR  KRAS  EGFR - KRAS        57       154 0.11354582 0.3067729
#> STK11 - TP53 STK11  TP53 STK11 - TP53        59       221 0.11752988 0.4402390
#> BRAF - KRAS   BRAF  KRAS  BRAF - KRAS        35       154 0.06972112 0.3067729
#> KRAS - STK11  KRAS STK11 KRAS - STK11       154        59 0.30677291 0.1175299
#> EGFR - TP53   EGFR  TP53  EGFR - TP53        57       221 0.11354582 0.4402390
#>              overlap  w_overlap max_overlap freq_overlap r_overlap w_r_overlap
#> KRAS - TP53       49 35.8760174         154   0.31818182    97.750    58.42521
#> EGFR - KRAS        0  0.0000000          57   0.00000000    29.500    15.70896
#> STK11 - TP53      13  9.5456386          59   0.22033898    38.250    22.31070
#> BRAF - KRAS        2  0.9821429          35   0.05714286    18.625    10.05346
#> KRAS - STK11      28 25.9230769          59   0.47457627    30.750    16.86851
#> EGFR - TP53       32 29.3335097          57   0.56140351    37.250    21.60725
#>                     wES wFDR        nES mean_r_nES nFDR cum_freq nFDR2 type
#> KRAS - TP53  -15.944687    0 -13.652182 -2.2925057    0      375     0   ME
#> EGFR - KRAS  -11.107910    0  -9.569699 -1.5382116    0      211     0   ME
#> STK11 - TP53  -9.026262    0  -7.230882 -1.7953795    0      280     0   ME
#> BRAF - KRAS   -6.414386    0  -5.265037 -1.1493495    0      189     0   ME
#> KRAS - STK11   6.402547    0   5.056438  1.3461097    0      213     0   CO
#> EGFR - TP53    5.463292    0   4.557487  0.9058047    0      278     0   CO
#>               FDR
#> KRAS - TP53  TRUE
#> EGFR - KRAS  TRUE
#> STK11 - TP53 TRUE
#> BRAF - KRAS  TRUE
#> KRAS - STK11 TRUE
#> EGFR - TP53  TRUE
# }
```
