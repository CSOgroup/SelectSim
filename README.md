
<!-- README.md is generated from README.Rmd. Please edit that file -->

# SelectSim <img src="man/figures/logo.png" align="right" height="138" alt="" />

<!-- badges: start -->
<!-- badges: end -->

The goal of SelectSim package is to implement the methodology to infer
inter-dependencies between functional alterations in cancer. SelectSim
estimates the expected number of mutations in a given gene and a given
sample from the mutation frequency of the gene, f(g), and the tumor
mutation burden (TMB) of the sample, $\mu$(t). These values can be
estimated within specific mutation and tumor subsets, to account for
heterogeneous tumor types, tissue specificities, and distinct mutational
processes

<figure>
<img src="vignettes/SelectSim_method.png" alt="SelectSim Method" />
<figcaption aria-hidden="true">SelectSim Method</figcaption>
</figure>

## Installation

- You can install the development version of SelectSim from
  [GitHub](https://github.com/CSOgroup/SelectSim) with:

``` r
# install.packages("devtools")
devtools::install_github("CSOgroup/SelectSim",dependencies = TRUE, build_vignettes = TRUE)
```

- For more details on installation refer to
  [INSTALLATION](INSTALLATION.md)

## Example

Check the the Vignettes folder or visit
[website](https://csogroup.github.io/SelectSim/)

### Who do I talk to?

- For any bugs or feature support in using SelectSim, please use the
  [issue tracker](https://github.com/CSOgroup/SelectSim/issues).
- For any other question related to SelectSim, please contact Prof
  Giovanni Ciriello (<giovanni.ciriello@unil.ch>) or Arvind Iyer
  (<arvind.iyer@unil.ch>).
