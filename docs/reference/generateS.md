# Generate S matrix

Generate S matrix

## Usage

``` r
generateS(gam, sample.weights, upperBound = 1)
```

## Arguments

- gam:

  the gam with genes\*samples

- sample.weights:

  the samples weights

- upperBound:

  clip the values greater than 1 to keep it bounded between 0 to 1

## Value

S the S matrix
