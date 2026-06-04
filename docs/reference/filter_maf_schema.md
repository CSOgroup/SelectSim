# This function filters a MAF dataframe by sample id

`filter_maf_schema()` takes a maf file and filters a MAF dataframe by
retaining only the rows with a column value included in the values list

## Usage

``` r
filter_maf_schema(maf, schema = TCGA_maf_schema, column, values, ...)
```

## Arguments

- maf:

  a maf as dataframe

- schema:

  a schema of datafrane check Select::TCGA_maf_schema for example

- column:

  column in maf file to filter

- values:

  a list containing the elements to file

- ...:

  Other options

## Value

filtered_maf a filtered maf file
