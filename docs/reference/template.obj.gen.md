# Generate the template matrix

Computes the expected mutation probability matrices (one per GAM type,
per sample block) used as the simulation background in
`null_model_parallel`.

## Usage

``` r
template.obj.gen(al)
```

## Arguments

- al:

  Alteration landscape object

## Value

A list with `template.obj` (per-block S matrices) and `temp_mat` (full
concatenated template matrices, one per GAM type).
