# Generate gam from the maf file

`maf2gam()` takes a maf file and converts into gam

## Usage

``` r
maf2gam(
  maf,
  sample.col = "Tumor_Sample_Barcode",
  gene.col = "Hugo_Symbol",
  value.var = "HGVSp_Short",
  samples = NULL,
  genes = NULL,
  fun.aggregate = length,
  binarize = TRUE,
  fill = NA,
  ...
)
```

## Arguments

- maf:

  A MAF dataframe

- sample.col:

  Column name for sample IDs

- gene.col:

  Column name for gene symbols

- value.var:

  Column used as the value to aggregate

- samples:

  Vector of sample IDs to include; NULL keeps all samples present in the
  MAF

- genes:

  Vector of gene names to include; NULL keeps all genes present in the
  MAF

- fun.aggregate:

  Aggregation function applied per (sample, gene) cell

- binarize:

  If TRUE, convert aggregated counts to binary presence/absence

- fill:

  Value used for missing (sample, gene) combinations

- ...:

  Additional arguments passed to reshape2::acast

## Value

Numeric matrix (samples x genes) representing the gene alteration
matrix.
