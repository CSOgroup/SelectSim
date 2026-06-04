# Filter a MAF dataframe by mutation type

Filter a MAF dataframe by mutation type

## Usage

``` r
filter_maf_mutation.type(
  maf,
  variants,
  variant.col = "Variant_Classification",
  ...
)
```

## Arguments

- maf:

  A MAF dataframe

- variants:

  Vector of variant classification values to retain

- variant.col:

  Column name containing variant classifications

- ...:

  Additional arguments passed to filter_maf_column

## Value

Filtered MAF dataframe.
