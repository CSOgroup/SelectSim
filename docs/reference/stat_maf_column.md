# Summary functions for MAF file

`stat_maf_column()` takes a maf file and filters a MAF dataframe by
retaining only the rows with a column value included in the values list

## Usage

``` r
stat_maf_column(maf, column, ...)
```

## Arguments

- maf:

  a maf as dataframe

- column:

  a schema of datafrane check Select::TCGA_maf_schema for example

- ...:

  Other options

## Value

filtered_maf a filtered maf file
