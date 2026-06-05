## R CMD check results

0 errors | 0 warnings | 0 notes

## Test environments

- local macOS (darwin 24.6.0), R 4.x
- GitHub Actions: ubuntu-latest (R release, R devel), macos-latest (R release), windows-latest (R release)

## Downstream dependencies

None on CRAN yet (new submission).

## Notes to CRAN reviewers

- This package includes compiled C++ code via Rcpp/RcppArmadillo.
- `tictoc` is in Suggests (not Imports); timing output is only shown when
  `verbose = TRUE` and the package is available.
- All examples that run the main `selectX()` function are wrapped in
  `\donttest{}` because the default permutation count (n.permut = 1000)
  takes ~30 s on a single core; the examples use n.permut = 10 to keep
  wall time short.
