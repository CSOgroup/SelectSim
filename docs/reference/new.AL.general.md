# Create an AL object

Create an Alteration Landscape (AL) object which contains gams and
mutation burden of samples of associated gams.

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

  The binary alteration matrix with row as features and column as
  samples.

- feat.covariates:

  gene/feature covariate.

- sample.covariates:

  sample covariate.

- min.freq:

  minimum frequency of genes to be mutated.

- verbose:

  print the time and each steps.

## Value

An Alteration Landscape (AL) object with the gam.
