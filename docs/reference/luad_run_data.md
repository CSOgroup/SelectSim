# Lung adenocarcinoma from TCGA cohort as SelectSim run object

Pre-processed TCGA LUAD data ready to pass directly to
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
