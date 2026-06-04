# Filter a MAF dataframe by specific gene-mutation combinations

Filter a MAF dataframe by specific gene-mutation combinations

## Usage

``` r
filter_maf_mutations(
  maf,
  values,
  maf.col = c("Hugo_Symbol", "HGVSp_Short"),
  values.col = maf.col,
  ...
)
```

## Arguments

- maf:

  A MAF dataframe

- values:

  Dataframe of allowed (gene, mutation) combinations

- maf.col:

  Columns in maf to join on

- values.col:

  Corresponding columns in values to join on

- ...:

  Additional arguments passed to filter_maf_complex

## Value

Filtered MAF dataframe containing only rows matching the allowed
combinations.
