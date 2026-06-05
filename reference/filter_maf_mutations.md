# Filter a MAF dataframe by specific gene-mutation combinations

Filter a MAF dataframe by specific gene-mutation combinations

## Usage

``` r
filter_maf_mutations(
  maf,
  values,
  maf.col = c("Hugo_Symbol", "HGVSp_Short"),
  values.col = maf.col,
  ...
)
```

## Arguments

- maf:

  A MAF dataframe

- values:

  Dataframe of allowed (gene, mutation) combinations

- maf.col:

  Columns in maf to join on

- values.col:

  Corresponding columns in values to join on

- ...:

  Additional arguments passed to filter_maf_complex

## Value

Filtered MAF dataframe containing only rows matching the allowed
combinations.

## Examples

``` r
data(luad_maf, package = "SelectSim")
allowed <- data.frame(Hugo_Symbol  = c("TP53", "KRAS"),
                      HGVSp_Short  = c("p.R175H", "p.G12C"))
filter_maf_mutations(luad_maf, allowed)
#>    Hugo_Symbol HGVSp_Short Chromosome Start_Position End_Position
#> 1         KRAS      p.G12C         12       25398285     25398285
#> 2         KRAS      p.G12C         12       25398285     25398285
#> 3         KRAS      p.G12C         12       25398285     25398285
#> 4         KRAS      p.G12C         12       25398285     25398285
#> 5         KRAS      p.G12C         12       25398285     25398285
#> 6         KRAS      p.G12C         12       25398285     25398285
#> 7         KRAS      p.G12C         12       25398285     25398285
#> 8         KRAS      p.G12C         12       25398285     25398285
#> 9         KRAS      p.G12C         12       25398285     25398285
#> 10        KRAS      p.G12C         12       25398285     25398285
#> 11        KRAS      p.G12C         12       25398285     25398285
#> 12        KRAS      p.G12C         12       25398285     25398285
#> 13        KRAS      p.G12C         12       25398285     25398285
#> 14        KRAS      p.G12C         12       25398285     25398285
#> 15        KRAS      p.G12C         12       25398285     25398285
#> 16        KRAS      p.G12C         12       25398285     25398285
#> 17        KRAS      p.G12C         12       25398285     25398285
#> 18        KRAS      p.G12C         12       25398285     25398285
#> 19        KRAS      p.G12C         12       25398285     25398285
#> 20        KRAS      p.G12C         12       25398285     25398285
#> 21        KRAS      p.G12C         12       25398285     25398285
#> 22        KRAS      p.G12C         12       25398285     25398285
#> 23        KRAS      p.G12C         12       25398285     25398285
#> 24        KRAS      p.G12C         12       25398285     25398285
#> 25        KRAS      p.G12C         12       25398285     25398285
#> 27        KRAS      p.G12C         12       25398285     25398285
#> 28        KRAS      p.G12C         12       25398285     25398285
#> 30        KRAS      p.G12C         12       25398285     25398285
#> 32        KRAS      p.G12C         12       25398285     25398285
#> 33        KRAS      p.G12C         12       25398285     25398285
#> 34        KRAS      p.G12C         12       25398285     25398285
#> 35        KRAS      p.G12C         12       25398285     25398285
#> 36        KRAS      p.G12C         12       25398285     25398285
#> 37        KRAS      p.G12C         12       25398285     25398285
#> 38        KRAS      p.G12C         12       25398285     25398285
#> 39        KRAS      p.G12C         12       25398285     25398285
#> 40        KRAS      p.G12C         12       25398285     25398285
#> 41        KRAS      p.G12C         12       25398285     25398285
#> 42        KRAS      p.G12C         12       25398285     25398285
#> 43        KRAS      p.G12C         12       25398285     25398285
#> 44        KRAS      p.G12C         12       25398285     25398285
#> 45        KRAS      p.G12C         12       25398285     25398285
#> 46        KRAS      p.G12C         12       25398285     25398285
#> 47        KRAS      p.G12C         12       25398285     25398285
#> 48        KRAS      p.G12C         12       25398285     25398285
#> 50        KRAS      p.G12C         12       25398285     25398285
#> 51        KRAS      p.G12C         12       25398285     25398285
#> 52        KRAS      p.G12C         12       25398285     25398285
#> 53        KRAS      p.G12C         12       25398285     25398285
#> 54        KRAS      p.G12C         12       25398285     25398285
#> 55        KRAS      p.G12C         12       25398285     25398285
#> 56        KRAS      p.G12C         12       25398285     25398285
#> 57        KRAS      p.G12C         12       25398285     25398285
#> 58        KRAS      p.G12C         12       25398285     25398285
#> 59        KRAS      p.G12C         12       25398285     25398285
#> 60        KRAS      p.G12C         12       25398285     25398285
#> 61        KRAS      p.G12C         12       25398285     25398285
#> 62        KRAS      p.G12C         12       25398285     25398285
#> 63        KRAS      p.G12C         12       25398285     25398285
#> 64        KRAS      p.G12C         12       25398285     25398285
#> 65        KRAS      p.G12C         12       25398285     25398285
#> 66        KRAS      p.G12C         12       25398285     25398285
#> 67        KRAS      p.G12C         12       25398285     25398285
#> 68        KRAS      p.G12C         12       25398285     25398285
#> 69        KRAS      p.G12C         12       25398285     25398285
#> 70        TP53     p.R175H         17        7578406      7578406
#> 71        TP53     p.R175H         17        7578406      7578406
#> 72        TP53     p.R175H         17        7578406      7578406
#>    Variant_Classification         Tumor_Sample_Barcode          sample
#> 1       Missense_Mutation TCGA-73-7498-01A-12D-2184-08 TCGA-73-7498-01
#> 2       Missense_Mutation TCGA-55-8299-01A-11D-2284-08 TCGA-55-8299-01
#> 3       Missense_Mutation TCGA-55-A490-01A-11D-A24D-08 TCGA-55-A490-01
#> 4       Missense_Mutation TCGA-55-8302-01A-11D-2323-08 TCGA-55-8302-01
#> 5       Missense_Mutation TCGA-64-5774-01A-01D-1625-08 TCGA-64-5774-01
#> 6       Missense_Mutation TCGA-78-7145-01A-11D-2036-08 TCGA-78-7145-01
#> 7       Missense_Mutation TCGA-44-7671-01A-11D-2063-08 TCGA-44-7671-01
#> 8       Missense_Mutation TCGA-55-1595-01A-01D-0969-08 TCGA-55-1595-01
#> 9       Missense_Mutation TCGA-75-5126-01A-01D-1753-08 TCGA-75-5126-01
#> 10      Missense_Mutation TCGA-93-7347-01A-11D-2184-08 TCGA-93-7347-01
#> 11      Missense_Mutation TCGA-86-7713-01A-11D-2063-08 TCGA-86-7713-01
#> 12      Missense_Mutation TCGA-55-8097-01A-11D-2238-08 TCGA-55-8097-01
#> 13      Missense_Mutation TCGA-05-4250-01A-01D-1105-08 TCGA-05-4250-01
#> 14      Missense_Mutation TCGA-05-4244-01A-01D-1105-08 TCGA-05-4244-01
#> 15      Missense_Mutation TCGA-91-6849-01A-11D-1945-08 TCGA-91-6849-01
#> 16      Missense_Mutation TCGA-97-A4M5-01A-11D-A24P-08 TCGA-97-A4M5-01
#> 17      Missense_Mutation TCGA-73-4662-01A-01D-1265-08 TCGA-73-4662-01
#> 18      Missense_Mutation TCGA-MP-A4TD-01A-32D-A25L-08 TCGA-MP-A4TD-01
#> 19      Missense_Mutation TCGA-05-4403-01A-01D-1265-08 TCGA-05-4403-01
#> 20      Missense_Mutation TCGA-05-4415-01A-22D-1855-08 TCGA-05-4415-01
#> 21      Missense_Mutation TCGA-69-7973-01A-11D-2184-08 TCGA-69-7973-01
#> 22      Missense_Mutation TCGA-95-7562-01A-11D-2238-08 TCGA-95-7562-01
#> 23      Missense_Mutation TCGA-55-7907-01A-11D-2167-08 TCGA-55-7907-01
#> 24      Missense_Mutation TCGA-4B-A93V-01A-11D-A397-08 TCGA-4B-A93V-01
#> 25      Missense_Mutation TCGA-55-6983-01A-11D-1945-08 TCGA-55-6983-01
#> 27      Missense_Mutation TCGA-MP-A4TE-01A-22D-A25L-08 TCGA-MP-A4TE-01
#> 28      Missense_Mutation TCGA-55-7281-01A-11D-2036-08 TCGA-55-7281-01
#> 30      Missense_Mutation TCGA-50-5936-01A-11D-1625-08 TCGA-50-5936-01
#> 32      Missense_Mutation TCGA-50-5933-01A-11D-1753-08 TCGA-50-5933-01
#> 33      Missense_Mutation TCGA-05-4405-01A-21D-1855-08 TCGA-05-4405-01
#> 34      Missense_Mutation TCGA-35-3615-01A-01D-1040-01 TCGA-35-3615-01
#> 35      Missense_Mutation TCGA-NJ-A4YI-01A-11D-A25L-08 TCGA-NJ-A4YI-01
#> 36      Missense_Mutation TCGA-05-4418-01A-01D-1265-08 TCGA-05-4418-01
#> 37      Missense_Mutation TCGA-64-1677-01A-01W-0928-08 TCGA-64-1677-01
#> 38      Missense_Mutation TCGA-55-7914-01A-11D-2167-08 TCGA-55-7914-01
#> 39      Missense_Mutation TCGA-MP-A4TI-01A-21D-A24P-08 TCGA-MP-A4TI-01
#> 40      Missense_Mutation TCGA-78-7539-01A-11D-2063-08 TCGA-78-7539-01
#> 41      Missense_Mutation TCGA-50-7109-01A-11D-2036-08 TCGA-50-7109-01
#> 42      Missense_Mutation TCGA-50-5932-01A-11D-1753-08 TCGA-50-5932-01
#> 43      Missense_Mutation TCGA-97-7938-01A-11D-2167-08 TCGA-97-7938-01
#> 44      Missense_Mutation TCGA-55-8207-01A-11D-2238-08 TCGA-55-8207-01
#> 45      Missense_Mutation TCGA-50-8459-01A-11D-2323-08 TCGA-50-8459-01
#> 46      Missense_Mutation TCGA-64-5778-01A-01D-1625-08 TCGA-64-5778-01
#> 47      Missense_Mutation TCGA-55-8203-01A-11D-2238-08 TCGA-55-8203-01
#> 48      Missense_Mutation TCGA-05-4417-01A-22D-1855-08 TCGA-05-4417-01
#> 50      Missense_Mutation TCGA-78-7148-01A-11D-2036-08 TCGA-78-7148-01
#> 51      Missense_Mutation TCGA-67-3774-01A-01D-1040-01 TCGA-67-3774-01
#> 52      Missense_Mutation TCGA-78-7167-01A-11D-2063-08 TCGA-78-7167-01
#> 53      Missense_Mutation TCGA-78-7166-01A-12D-2063-08 TCGA-78-7166-01
#> 54      Missense_Mutation TCGA-50-5051-01A-21D-1855-08 TCGA-50-5051-01
#> 55      Missense_Mutation TCGA-99-8032-01A-11D-2238-08 TCGA-99-8032-01
#> 56      Missense_Mutation TCGA-44-7659-01A-11D-2063-08 TCGA-44-7659-01
#> 57      Missense_Mutation TCGA-55-7726-01A-11D-2167-08 TCGA-55-7726-01
#> 58      Missense_Mutation TCGA-MP-A4TF-01A-11D-A25L-08 TCGA-MP-A4TF-01
#> 59      Missense_Mutation TCGA-05-4249-01A-01D-1105-08 TCGA-05-4249-01
#> 60      Missense_Mutation TCGA-55-8508-01A-11D-2393-08 TCGA-55-8508-01
#> 61      Missense_Mutation TCGA-97-A4M0-01A-11D-A24P-08 TCGA-97-A4M0-01
#> 62      Missense_Mutation TCGA-NJ-A55O-01A-11D-A25L-08 TCGA-NJ-A55O-01
#> 63      Missense_Mutation TCGA-55-8616-01A-11D-2393-08 TCGA-55-8616-01
#> 64      Missense_Mutation TCGA-49-4505-01A-01D-1931-08 TCGA-49-4505-01
#> 65      Missense_Mutation TCGA-55-8615-01A-11D-2393-08 TCGA-55-8615-01
#> 66      Missense_Mutation TCGA-J2-8194-01A-11D-2238-08 TCGA-J2-8194-01
#> 67      Missense_Mutation TCGA-99-8028-01A-11D-2238-08 TCGA-99-8028-01
#> 68      Missense_Mutation TCGA-MN-A4N5-01A-11D-A24P-08 TCGA-MN-A4N5-01
#> 69      Missense_Mutation TCGA-69-7978-01A-11D-2184-08 TCGA-69-7978-01
#> 70      Missense_Mutation TCGA-49-4514-01A-21D-1855-08 TCGA-49-4514-01
#> 71      Missense_Mutation TCGA-05-4405-01A-21D-1855-08 TCGA-05-4405-01
#> 72      Missense_Mutation TCGA-50-6595-01A-12D-1855-08 TCGA-50-6595-01
```
