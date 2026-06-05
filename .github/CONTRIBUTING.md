# Contributing to SelectSim

Thank you for your interest in contributing to SelectSim. This document explains how to report bugs, suggest features, and submit code.

## Reporting bugs

Please use the [issue tracker](https://github.com/CSOgroup/SelectSim/issues) and select the **Bug report** template. Include:

- SelectSim version (`packageVersion("SelectSim")`)
- R version and OS
- A minimal reproducible example (preferably using the bundled `luad_run_data`)
- The full error message or unexpected output

## Suggesting features

Open an issue using the **Feature request** template. Describe the use case and why it would benefit users of the package.

## Contributing code

1. Fork the repository and create a branch from `main`.
2. Follow the existing code style (base R where possible; tidyverse only via already-imported packages).
3. Add or update tests in `tests/testthat/` for any new or changed behaviour.
4. Run `devtools::check()` and ensure 0 errors, 0 warnings, and 0 notes before opening a pull request.
5. Update roxygen2 documentation (`devtools::document()`) if you change function signatures or add new functions.
6. Open a pull request against `main` using the provided PR template.

## Code of conduct

This project follows the [Contributor Covenant Code of Conduct](CODE_OF_CONDUCT.md). By participating you agree to abide by its terms.

## Questions

For questions about the SelectSim methodology, contact Giovanni Ciriello (giovanni.ciriello@unil.ch) or Arvind Iyer (ayalurarvind@gmail.com).
