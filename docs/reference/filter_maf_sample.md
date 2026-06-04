# Filter a MAF dataframe by sample ID

Filter a MAF dataframe by sample ID

## Usage

``` r
filter_maf_sample(maf, samples, sample.col = "Tumor_Sample_Barcode", ...)
```

## Arguments

- maf:

  A MAF dataframe

- samples:

  Vector of sample IDs to retain

- sample.col:

  Column name containing sample IDs

- ...:

  Additional arguments passed to filter_maf_column

## Value

Filtered MAF dataframe.
