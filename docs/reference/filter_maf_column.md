# Filter maf function

`filter_maf_column()` takes a maf file and filters a MAF dataframe by
retaining only the rows with a column value included in the values list

## Usage

``` r
filter_maf_column(maf, values, column, inclusive = TRUE, fixed = TRUE, ...)
```

## Arguments

- maf:

  a maf as dataframe

- values:

  a list containing the elements to filter

- column:

  column in maf file to filter

- inclusive:

  a boolena to include or exclude the dataframe with values in list
  provided

- fixed:

  a grep argument to specify if grep use the argumnet as string or not

- ...:

  Other options

## Value

filtered_maf a filtered maf file
