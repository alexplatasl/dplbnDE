## Resubmission

This package was previously archived from CRAN. This version (0.2.0) addresses
the issues that led to its removal and includes several improvements:

- Fixed compatibility with current bnclassify (>= 0.4.5) on CRAN
- Fixed bug in structure selection that ignored user-specified structures
- Fixed custom Bayesian network input via edgelist parameter
- Performance optimizations in the inner DE loop
- Added matrixStats as dependency

## Test environments

* Windows 11 (local), R 4.5.0
* R CMD check results: 0 errors | 0 warnings | 1 note

## R CMD check results

There was 1 NOTE:

* checking for future file timestamps ... NOTE
  unable to verify current time

This is a transient network issue unrelated to the package.

## Downstream dependencies

There are currently no downstream dependencies for this package.
