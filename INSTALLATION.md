# Installation

## Prerequisites

SelectSim requires **R ≥ 3.5** and includes compiled C++ code via `Rcpp` and `RcppArmadillo`.
Before installing, make sure you have a working C++ compiler:

- **macOS** — Install Xcode Command Line Tools: `xcode-select --install`
- **Linux** — Install build tools, e.g. on Ubuntu/Debian: `sudo apt install build-essential`
- **Windows** — Install [Rtools](https://cran.r-project.org/bin/windows/Rtools/) and ensure it is on your PATH

---

## Option 1 — Install from GitHub (standard)

Install directly from GitHub using `pak`:

```r
# install.packages("pak")
pak::pak("CSOgroup/SelectSim")
```

Or using `devtools` (legacy):

```r
# install.packages("devtools")
devtools::install_github("CSOgroup/SelectSim", dependencies = TRUE, build_vignettes = TRUE)
```

---

## Option 2 — Install with micromamba (recommended)

[micromamba](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html) provides a fully isolated, reproducible R environment. This is the **recommended approach** if you are on a shared server, HPC cluster, or want to avoid conflicts with your existing R installation.

### Step 1 — Install micromamba

Follow the [official micromamba installation guide](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html) for your platform.

### Step 2 — Create and activate an R environment

```bash
micromamba create -n r_env
micromamba activate r_env
```

### Step 3 — Install R and required system dependencies

```bash
micromamba install \
  conda-forge::r-base \
  conda-forge::r-essentials \
  conda-forge::r-devtools \
  conda-forge::r-pak \
  conda-forge::cmake \
  conda-forge::r-rcppparallel \
  conda-forge::r-rfast
```

### Step 4 — Install SelectSim inside R

```r
# install.packages("pak")
pak::pak("CSOgroup/SelectSim")
```

---

## Option 3 — Restore the exact development environment

A fully pinned conda environment file is provided at [`dev/environment.yml`](dev/environment.yml).
This reproduces the exact software stack (R version, all system libraries) used during development.

```bash
micromamba create -f dev/environment.yml
micromamba activate r_env
```

Then install SelectSim:

```r
pak::pak("CSOgroup/SelectSim")
```

> **Note:** The environment file pins exact package versions and is Linux-specific (linux-64).
> It may not work on macOS or Windows.

---

## Verifying the installation

After installing, confirm everything works with a quick test run:

```r
library(SelectSim)
data(luad_run_data, package = "SelectSim")

result <- selectX(
  M                = luad_run_data$M,
  sample.class     = luad_run_data$sample.class,
  alteration.class = luad_run_data$alteration.class,
  n.cores          = 1,
  min.freq         = 10,
  n.permut         = 10    # small number for a quick smoke test
)
```

A successful run returns a list with a `result` data frame of evolutionary dependencies.
For a full analysis, use `n.permut = 1000` or higher.

---

## Troubleshooting

**`rfast` fails to install**
`rfast` requires a Fortran compiler.
On macOS, install `gfortran` from [mac.r-project.org/tools](https://mac.r-project.org/tools/).
On Linux, install via your package manager (e.g. `sudo apt install gfortran`).
With micromamba, `conda-forge::r-rfast` handles this automatically.

**Vignettes not building**
Ensure `pandoc` is installed (`pandoc --version` in a terminal).
Install via `brew install pandoc` (macOS), `sudo apt install pandoc` (Linux), or through the micromamba environment (`conda-forge::pandoc`).

**Compilation errors on Windows**
Ensure Rtools is installed and its `bin/` directory is on your PATH.
Run `pkgbuild::check_build_tools()` in R to verify your toolchain is detected correctly.

**Vignettes missing after install**
`pak::pak()` does not build vignettes by default. If you need the vignettes locally, use:
```r
devtools::install_github("CSOgroup/SelectSim", dependencies = TRUE, build_vignettes = TRUE)
```
