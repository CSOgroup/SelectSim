# Generating the null_simulation matrix

Generating the null_simulation matrix

## Usage

``` r
null_model_parallel_debug(al, temp_mat, W, n.cores = 1, n.permut)
```

## Arguments

- al:

  Alteration landscape object

- temp_mat:

  template matrixes

- W:

  weight matrix

- n.cores:

  Number of cores

- n.permut:

  Number of simulation

## Value

Template matrix as list object

Ordering of genes impacts the results as residual subraction is not
correct.
