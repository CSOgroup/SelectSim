# Extract null-model weighted overlap distribution for a gene pair

Returns the vector of weighted co-mutation values from all null-model
permutations for a specific pair of genes.

## Usage

``` r
overlap_pair_extract(gene1, gene2, obj)
```

## Arguments

- gene1:

  Name of the first gene/alteration.

- gene2:

  Name of the second gene/alteration.

- obj:

  SelectX object (list returned by `selectX()$obj`).

## Value

Numeric vector of length `obj$nSim` with the null-model weighted overlap
values for the pair.

## Examples

``` r
# \donttest{
data(luad_run_data, package = "SelectSim")
result <- selectX(M = luad_run_data$M,
                  sample.class = luad_run_data$sample.class,
                  alteration.class = luad_run_data$alteration.class,
                  n.cores = 1, min.freq = 10, n.permut = 10,
                  verbose = FALSE)
genes <- rownames(result$obj$al$am$full)[1:2]
overlap_pair_extract(genes[1], genes[2], result$obj)
#> [1] 0.7692308 1.3942308 1.6526894 0.9821429 0.6451613 1.1692308 1.0918114
#> [8] 0.4000000
# }
```
