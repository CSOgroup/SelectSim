# Create an Alteration Landscape (AL) object

Builds an Alteration Landscape object from a list of genome alteration
matrices and their corresponding tumor mutation burdens.

## Usage

``` r
new.AL.general(
  am,
  feat.covariates = NULL,
  sample.covariates = NULL,
  min.freq,
  verbose = FALSE
)
```

## Arguments

- am:

  A named list with two required elements: `M` (a named list of binary
  alteration matrices, each genes x samples) and `tmb` (a named list of
  data frames, one per matrix in `M`, each with columns `sample` and
  `mutation`). The names of `M` and `tmb` must match. All matrices in
  `M` must have identical row and column names.

- feat.covariates:

  Named character vector of alteration-type annotations, one entry per
  feature (gene). Names must match rownames of the matrices in `M`. If
  `NULL`, all features are labelled `"MUT"`.

- sample.covariates:

  Named character vector of sample-type annotations, one entry per
  sample. Names must match colnames of the matrices in `M`. If `NULL`,
  all samples are labelled `"sample"`.

- min.freq:

  Minimum number of samples a gene must be mutated in (strictly greater
  than) to be retained. Features with `rowSums <= min.freq` are dropped.

- verbose:

  Logical; print progress messages.

## Value

An Alteration Landscape (AL) object (list of class `"AL"`) containing
the filtered alteration matrices, TMB vectors, and covariate
assignments.

## Examples

``` r
# \donttest{
data(luad_run_data, package = "SelectSim")
al <- new.AL.general(
  am               = luad_run_data$M,
  feat.covariates  = luad_run_data$alteration.class,
  sample.covariates = luad_run_data$sample.class,
  min.freq         = 10
)
# }
```
