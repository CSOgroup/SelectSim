# Generating the weight matrix

Generating the weight matrix

## Usage

``` r
generateW_mean_tmb(
  tmb,
  mean_tmb,
  ngenes,
  lambda = 0.3,
  tau = 1,
  discrete = TRUE
)
```

## Arguments

- tmb:

  TMB dataframe

- mean_tmb:

  TMB dataframe

- ngenes:

  Number of genes

- lambda:

  weight penalty factor (default 0.3)

- tau:

  fold change factor

- discrete:

  True discrete weights

## Value

weight matrix (i.e vector as matrix)
