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

## Examples

``` r
# \donttest{
data(luad_run_data, package = "SelectSim")
al       <- new.AL.general(luad_run_data$M,
                           feat.covariates  = luad_run_data$alteration.class,
                           sample.covariates = luad_run_data$sample.class,
                           min.freq = 10)
temp_obj <- template.obj.gen(al)
W        <- generateW_block(al, lambda = 0.3, tau = 1)
sims     <- null_model_parallel(al, temp_obj$temp_mat, W$W,
                                n.cores = 1, n.permut = 10)
length(sims)
#> [1] 10
# }
```
