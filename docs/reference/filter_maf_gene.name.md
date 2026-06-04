# Filter a MAF dataframe by gene name

Filter a MAF dataframe by gene name

## Usage

``` r
filter_maf_gene.name(maf, genes, gene.col = "Hugo_Symbol", ...)
```

## Arguments

- maf:

  A MAF dataframe

- genes:

  Vector of gene names to retain

- gene.col:

  Column name containing gene symbols

- ...:

  Additional arguments passed to filter_maf_column

## Value

Filtered MAF dataframe.
