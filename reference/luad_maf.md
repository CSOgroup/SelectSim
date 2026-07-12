# Lung adenocarcinoma MAF from TCGA cohort

MAF file of LUAD from TCGA

## Usage

``` r
data(luad_maf)
```

## Format

A data frame

## Value

A data frame containing TCGA LUAD mutation data.

## Examples

``` r
data(luad_maf)
dim(luad_maf)
#> [1] 220734      8
head(luad_maf)
#>   Chromosome Start_Position End_Position   Hugo_Symbol Variant_Classification
#> 1         10      101814119    101814119          CPN1      Missense_Mutation
#> 2         10      129902901    129902901         MKI67                 Silent
#> 3         10       21104601     21104606          NEBL           In_Frame_Del
#> 4         10       45652518     45652518 RP11-445N18.7                    RNA
#> 5         10       50667200     50667200         ERCC6                 Silent
#> 6         10         532472       532472         DIP2C            Splice_Site
#>           Tumor_Sample_Barcode          sample    HGVSp_Short
#> 1 TCGA-05-4244-01A-01D-1105-08 TCGA-05-4244-01        p.H366D
#> 2 TCGA-05-4244-01A-01D-1105-08 TCGA-05-4244-01       p.N2401N
#> 3 TCGA-05-4244-01A-01D-1105-08 TCGA-05-4244-01 p.S730_V731del
#> 4 TCGA-05-4244-01A-01D-1105-08 TCGA-05-4244-01              .
#> 5 TCGA-05-4244-01A-01D-1105-08 TCGA-05-4244-01       p.S1381S
#> 6 TCGA-05-4244-01A-01D-1105-08 TCGA-05-4244-01   p.X29_splice
```
