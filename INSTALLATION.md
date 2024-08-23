## Installation

- You can install the development version of SelectSim from
  [GitHub](https://github.com/CSOgroup/SelectSim) with:

``` r
# install.packages("devtools")
devtools::install_github("CSOgroup/SelectSim",dependencies = TRUE, build_vignettes = TRUE)
```

## Installation with micromamba enviorment


`micromamba` is a tiny version of the mamba package manager (Like Conda). Refer [website](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html) for it's installation guide.

Steps to follow after installing micromamba in a terminal:

`micromamba create -n r_env r-essentials r-base`
`micromamba activate r_env`
`micromamba install conda-forge::r-devtools`
`micromamba install conda-forge::armadillo`
`micromamba install conda-forge::r-rcppparallel`
`micromamba install conda-forge::r-rfast`

- After this run R and install `SelectSim` R package as follows
``` r
# install.packages("devtools")
devtools::install_github("CSOgroup/SelectSim",dependencies = TRUE, build_vignettes = TRUE)
```

Alternative way install with provided `enviorment.yml` file.

`micromamba create -f env.yml`
`micromamba activate r_env`
- After this run R and install `SelectSim` R package as follows
``` r
# install.packages("devtools")
devtools::install_github("CSOgroup/SelectSim",dependencies = TRUE, build_vignettes = TRUE)
```

