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

## Examples

``` r
data(luad_maf, package = "SelectSim")
head(stat_maf_gene(luad_maf))
#> v
#>    A1BG    A1CF     A2M   A2ML1 A3GALT2  A4GALT 
#>      17      20      26      36       1       1 
```
