# Generate a pairs background plot for two dataset comaprision

Create an Alteration Landscape (AL) object which contains gams and
mutation burden of samples of associated gams.

## Usage

``` r
ridge_plot_ed_compare(result_df, obj1, obj2, name1, name2)
```

## Arguments

- result_df:

  a common result table of selectX obj\$result with dataset1_w_overlap
  and dataset2_w_overlap columns

- obj1:

  selectX obj of dataset1

- obj2:

  selectX obj of dataset2

- name1:

  dataset1 name

- name2:

  dataset2 name

## Value

a ggplot2 object with observed vs random overlap plot.
