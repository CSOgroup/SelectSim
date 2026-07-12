# Lung adenocarcinoma from TCGA cohort as SelectSim run object

Preprocessed TCGA LUAD data ready to pass directly to
[`selectX()`](https://csogroup.github.io/SelectSim/reference/selectX.md).

## Usage

``` r
data(luad_run_data)
```

## Format

A named list with four elements:

- M:

  A named list containing:

  M

  :   A named list of binary alteration matrices (genes x samples), one
      per alteration type (e.g., `missense`, `truncating`).

  tmb

  :   A named list of data frames, one per alteration type, each with
      columns `sample` (character) and `mutation` (integer TMB count).

- sample.class:

  Named character vector of sample-type annotations (length = number of
  samples). Names are sample IDs.

- alteration.class:

  Named character vector of alteration-type annotations (length = number
  of genes). Names are gene symbols.

## Value

A named list containing preprocessed LUAD input data for
[`selectX()`](https://csogroup.github.io/SelectSim/reference/selectX.md).

## Examples

``` r
data(luad_run_data)
names(luad_run_data)
#> [1] "M"                "sample.class"     "alteration.class"
str(luad_run_data, max.level = 1)
#> List of 3
#>  $ M               :List of 2
#>  $ sample.class    : Named chr [1:502] "LUAD" "LUAD" "LUAD" "LUAD" ...
#>   ..- attr(*, "names")= chr [1:502] "TCGA-05-4244-01" "TCGA-05-4249-01" "TCGA-05-4250-01" "TCGA-05-4382-01" ...
#>  $ alteration.class: Named chr [1:396] "MUT" "MUT" "MUT" "MUT" ...
#>   ..- attr(*, "names")= chr [1:396] "AKT1" "ALK" "APC" "AR" ...
```
