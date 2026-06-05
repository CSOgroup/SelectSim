# Changelog

## SelectSim 0.1.6

- Removed `reshape2` dependency; replaced with base R equivalents
  (`tapply`, `as.data.frame.table`).
- Fixed `theme_ridges()` partial argument match (`center` →
  `center_axis_labels`).
- Fixed parallel log file leak in `null_model_parallel` using
  [`tempfile()`](https://rdrr.io/r/base/tempfile.html) and
  [`on.exit()`](https://rdrr.io/r/base/on.exit.html).
- Refactored
  [`selectX()`](https://csogroup.github.io/SelectSim/reference/selectX.md)
  to eliminate duplicate verbose/non-verbose code paths.
- Added
  [`utils::globalVariables()`](https://rdrr.io/r/utils/globalVariables.html)
  declarations to suppress R CMD check NOTEs for ggplot2 aesthetics.
- Dropped `SystemRequirements: C++14` from DESCRIPTION.
- Added GitHub community health files: `CONTRIBUTING.md`,
  `CODE_OF_CONDUCT.md`, issue templates, PR template, and CI workflow.
- Expanded test coverage with `test-gam_utils.R`.
- Repository migrated: development moved from personal branches to
  [CSOgroup/SelectSim](https://github.com/CSOgroup/SelectSim); `dev`
  branch merged into `main` for initial public release.
- Bumped version to 0.1.6 across all source files; updated `DESCRIPTION`
  date to 2026-06-05.

## SelectSim 0.0.1.3

- Added a `NEWS.md` file to track changes to the package.
- Rename the pacakge SelectSim to SelectSim by creating a new git
  folder.
- Create CSO group repositroy and move the code there.
- Created the website for github

## SelectSim 0.0.1.4

- Added Mijan in author’s list
- Remove C/C++ code dependecny to avoid installation diffculties in
  different systems.
  - Hence move to using `Matrix` library functions and removed `RCpp`
    functions and code.
- Update the website and vignette accordingly.

## SelectSim 0.0.1.5

- Fixed bug of outlier functions and added C/C++ code back.
- More description in functions
