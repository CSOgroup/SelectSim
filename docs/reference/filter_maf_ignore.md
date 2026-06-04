# This function filters a MAF dataframe by retaining (or discarding) ignore mutations

`filter_maf_ignore()` takes a maf file and filters a MAF dataframe by
retaining only the rows with a column value included in the values list

## Usage

``` r
filter_maf_ignore(maf, schema = TCGA_maf_schema, ...)
```

## Arguments

- maf:

  a maf as dataframe

- schema:

  a schema of datafrane check Select::TCGA_maf_schema for example

- ...:

  Other options

## Value

filtered_maf a filtered maf file
