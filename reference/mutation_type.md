# Mutation list object

Mutation list object

## Usage

``` r
mutation_type
```

## Value

A list containing the supported mutation-type classifications.

## Examples

``` r
str(mutation_type)
#> List of 3
#>  $ truncating: chr [1:6] "Nonsense_Mutation" "Frame_Shift_Ins" "Frame_Shift_Del" "Splice_Site" ...
#>  $ missense  : chr [1:2] "Missense_Mutation" "Splice_Site"
#>  $ ignore    : chr [1:16] "Silent" "lincRNA" "IGR" "3'UTR" ...
```
