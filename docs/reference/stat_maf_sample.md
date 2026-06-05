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

## Examples

``` r
data(luad_maf, package = "SelectSim")
head(stat_maf_sample(luad_maf))
#> v
#> TCGA-05-4244-01A-01D-1105-08 TCGA-05-4249-01A-01D-1105-08 
#>                          270                          413 
#> TCGA-05-4250-01A-01D-1105-08 TCGA-05-4382-01A-01D-1931-08 
#>                          461                         2155 
#> TCGA-05-4384-01A-01D-1753-08 TCGA-05-4389-01A-01D-1265-08 
#>                          174                          305 
```
