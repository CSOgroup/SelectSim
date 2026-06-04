# Count mutations per sample in a MAF file

Count mutations per sample in a MAF file

## Usage

``` r
stat_maf_sample(maf, column = "Tumor_Sample_Barcode", ...)
```

## Arguments

- maf:

  A MAF dataframe

- column:

  Column name containing sample IDs

- ...:

  Additional arguments passed to stat_maf_column

## Value

Table of mutation counts per sample.
