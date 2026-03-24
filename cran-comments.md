## Resubmission

This is a resubmission addressing the issues from the previous review:

- Fixed typo "mutarion" -> "mutation" in DESCRIPTION
- Reduced example sizes (NP and G) so all examples run under 5s
- The words flagged as misspelled (Discriminative, Fukunaga, Sanderson,
  Storn, Tanabe, Zhang, independencies) are proper nouns, author names,
  and technical terms that are correctly spelled.

This package was previously archived from CRAN because it depended on
bnclassify which was also archived. bnclassify is back on CRAN (0.4.8).

## Test environments

* Windows 11 (local), R 4.5.0
* win-builder (r-devel-windows-x86_64) - previous submission
* Debian (r-devel-linux-x86_64-debian-gcc) - previous submission

## R CMD check results

0 errors | 0 warnings | 1 note

The NOTE about "New submission / Package was archived on CRAN" is expected
since this package is being resubmitted after archival.

## Downstream dependencies

There are currently no downstream dependencies for this package.
