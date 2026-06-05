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

## Examples

``` r
data(luad_maf, package = "SelectSim")
stat_maf_column(luad_maf, column = "Variant_Classification")
#> v
#>                3'Flank                  3'UTR                5'Flank 
#>                    624                   7860                    797 
#>                  5'UTR        Frame_Shift_Del        Frame_Shift_Ins 
#>                   4717                   3953                   1163 
#>           In_Frame_Del           In_Frame_Ins                 Intron 
#>                    373                     37                   5248 
#>      Missense_Mutation      Nonsense_Mutation       Nonstop_Mutation 
#>                 131221                  10905                    176 
#>                    RNA                 Silent            Splice_Site 
#>                   2980                  46060                   4399 
#> Translation_Start_Site 
#>                    221 
```
