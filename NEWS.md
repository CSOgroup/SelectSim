# SelectSim 0.1.6

* Removed `reshape2` dependency; replaced with base R equivalents (`tapply`, `as.data.frame.table`).
* Fixed `theme_ridges()` partial argument match (`center` → `center_axis_labels`).
* Fixed parallel log file leak in `null_model_parallel` using `tempfile()` and `on.exit()`.
* Refactored `selectX()` to eliminate duplicate verbose/non-verbose code paths.
* Added `utils::globalVariables()` declarations to suppress R CMD check NOTEs for ggplot2 aesthetics.
* Dropped `SystemRequirements: C++14` from DESCRIPTION.
* Added GitHub community health files: `CONTRIBUTING.md`, `CODE_OF_CONDUCT.md`, issue templates, PR template, and CI workflow.
* Expanded test coverage with `test-gam_utils.R`.
* Repository migrated: development moved from personal branches to [CSOgroup/SelectSim](https://github.com/CSOgroup/SelectSim); `dev` branch merged into `main` for initial public release.
* Bumped version to 0.1.6 across all source files
* Updated citation (`inst/CITATION`, `README`) to the published Nature Genetics article: Iyer A, Petrovic M, Sesia D, Nanni L, Mina M, Ciriello G (2026). Evolving patterns of co-mutations from tumor initiation to metastatic progression. *Nature Genetics*. DOI: 10.1038/s41588-026-02661-4

# SelectSim 0.0.1.3

* Added a `NEWS.md` file to track changes to the package.
* Renamed the package to SelectSim by creating a new Git repository.
* Created the CSOgroup repository and moved the code there.
* Created the GitHub website.

# SelectSim 0.0.1.4

* Added Mijan to the author list.
* Removed the C/C++ code dependency to avoid installation difficulties across systems.
  * Moved to `Matrix` package functions and removed the `Rcpp` functions and compiled code.
* Updated the website and vignette accordingly. 

# SelectSim 0.0.1.5

* Fixed bug of outlier functions and added C/C++ code back.
* More description in functions