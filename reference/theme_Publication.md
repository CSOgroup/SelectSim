# A clean ggplot2 theme for publication-quality plots

A clean ggplot2 theme for publication-quality plots

## Usage

``` r
theme_Publication(base_size = 14, base_family = "sans")
```

## Arguments

- base_size:

  size of fonts

- base_family:

  type of font

## Value

A ggplot2 theme object.

## Examples

``` r
library(ggplot2)
ggplot(data.frame(x = 1:3, y = 1:3), aes(x, y)) +
  geom_point() + theme_Publication()

```
