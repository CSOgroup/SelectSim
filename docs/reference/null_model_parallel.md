# Generating the null_simulation matrix

Generating the null_simulation matrix

## Usage

``` r
null_model_parallel(al, temp_mat, W, n.cores = 1, n.permut, seed = 42)
```

## Arguments

- al:

  Alteration landscape object

- temp_mat:

  template matrices

- W:

  weight matrix

- n.cores:

  Number of cores

- n.permut:

  Number of simulations

- seed:

  Random seed

## Value

List of simulated binary matrices (one per permutation)
