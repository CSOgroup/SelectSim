# SelectX main function to create alteration object with background model

[`selectX()`](https://csogroup.github.io/SelectSim/reference/selectX.md)
takes a list object which consist of genome alteration matrix and tumor
mutation burden.

## Usage

``` r
selectX_debug(
  M,
  sample.class,
  alteration.class,
  n.cores = 1,
  min.freq = 10,
  n.permut = 1000,
  lambda = 0.3,
  tao = 1,
  save.object = FALSE,
  folder = "./",
  verbose = TRUE,
  estimate_pairwise = FALSE,
  maxFDR = 0.25,
  seed = 42
)
```

## Arguments

- M:

  a list data object which consist of gams and tmbs.

- sample.class:

  sample covariates as named list.

- alteration.class:

  alteration covariates as named list.

- n.cores:

  no of cores.

- min.freq:

  number of samples for features to be atleast mutated in.

- n.permut:

  number of simulations.

- lambda:

  lambda parameter.

- tao:

  tao parameter.

- save.object:

  store the SelectX object.

- folder:

  folder path to store the results.

- verbose:

  print the time and each steps.

- estimate_pairwise:

  Compute pairwise p-value

- maxFDR:

  FDR value

- seed:

  a random seed

## Value

result a SelectX object with background model and other info along with
result table
