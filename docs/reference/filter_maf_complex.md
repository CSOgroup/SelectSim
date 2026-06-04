# Filter a MAF dataframe by a combination of column values

Filter a MAF dataframe by a combination of column values

## Usage

``` r
filter_maf_complex(maf, values, ...)
```

## Arguments

- maf:

  A MAF dataframe

- values:

  A dataframe of (column, value) pairs to match against

- ...:

  Additional arguments passed to merge

## Value

Filtered MAF dataframe containing only rows matching the value
combinations.
