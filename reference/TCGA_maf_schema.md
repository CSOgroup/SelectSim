# TCGA_maf_schema: schema for TCGA maf file to process the mutations

TCGA_maf_schema: schema for TCGA maf file to process the mutations

## Usage

``` r
TCGA_maf_schema
```

## Value

A list defining the expected TCGA MAF column schema.

## Examples

``` r
str(TCGA_maf_schema)
#> List of 3
#>  $ name         : chr "TCGA_maf"
#>  $ column       :List of 6
#>   ..$ gene         : chr "Hugo_Symbol"
#>   ..$ gene.name    : chr "Hugo_Symbol"
#>   ..$ sample       : chr "Tumor_Sample_Barcode"
#>   ..$ sample.name  : chr "Tumor_Sample_Barcode"
#>   ..$ mutation.type: chr "Variant_Classification"
#>   ..$ mutation     : chr "HGVSp_Short"
#>  $ mutation.type:List of 3
#>   ..$ truncating: chr [1:6] "Nonsense_Mutation" "Frame_Shift_Ins" "Frame_Shift_Del" "Splice_Site" ...
#>   ..$ missense  : chr [1:4] "Missense_Mutation" "Splice_Site" "In_Frame_Ins" "In_Frame_Del"
#>   ..$ ignore    : chr [1:16] "Silent" "lincRNA" "IGR" "3'UTR" ...
```
