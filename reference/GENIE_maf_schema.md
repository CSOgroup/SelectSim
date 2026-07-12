# GENIE_maf_schema: schema for GENIE maf file to process the mutations

GENIE_maf_schema: schema for GENIE maf file to process the mutations

## Usage

``` r
GENIE_maf_schema
```

## Value

A list defining the expected GENIE MAF column schema.

## Examples

``` r
str(GENIE_maf_schema)
#> List of 2
#>  $ column       :List of 6
#>   ..$ gene         : chr "Hugo_Symbol"
#>   ..$ gene.name    : chr "Hugo_Symbol"
#>   ..$ sample       : chr "Tumor_Sample_Barcode"
#>   ..$ sample.name  : chr "Tumor_Sample_Barcode"
#>   ..$ mutation.type: chr "Variant_Classification"
#>   ..$ mutation     : chr "HGVSp_Short"
#>  $ mutation.type:List of 3
#>   ..$ truncating: chr [1:6] "Nonsense_Mutation" "Frame_Shift_Ins" "Frame_Shift_Del" "Splice_Site" ...
#>   ..$ missense  : chr [1:2] "Missense_Mutation" "Splice_Site"
#>   ..$ ignore    : chr [1:16] "Silent" "lincRNA" "IGR" "3'UTR" ...
```
