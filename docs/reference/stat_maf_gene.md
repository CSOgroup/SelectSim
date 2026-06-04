# Count mutations per gene in a MAF file

Count mutations per gene in a MAF file

## Usage

``` r
stat_maf_gene(maf, column = "Hugo_Symbol", ...)
```

## Arguments

- maf:

  A MAF dataframe

- column:

  Column name containing gene symbols

- ...:

  Additional arguments passed to stat_maf_column

## Value

Table of mutation counts per gene.
