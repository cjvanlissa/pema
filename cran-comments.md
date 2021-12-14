# Version 0.1.0

* Addressed all comments by CRAN maintainer Julia Haider:
* Fixed comment by Julia Haider: Add references describing the methods 
  in the description field of your DESCRIPTION file
    + Added: 
  "see Van Lissa & Van Erp (2021). <doi:10.31234/osf.io/6phs5>."
* Fixed comment by Julia Haider: Add \value to .Rd files
    + Added \value to brma(), as.stan(), and maxap(). simulate_smd() already had
      a \value entry.

## Test environments
* local x86_64-pc-linux-gnu, R 4.1.2
* rhub::check_for_cran()
* devtools::check_release()
* devtools::check_oldrelease()
* devtools::check_devel()
* GitHub actions, windows-latest (release)
* GitHub actions, macOS-latest (release)
* GitHub actions, ubuntu-20.04 (release)
* GitHub actions, ubuntu-20.04 (devel)

## R CMD check results

0 errors | 0 warnings | 3 notes

* New submission
    + This is indeed a new submission
* GNU make is a SystemRequirements.
    + This package relies on rstan, which requires GNU make
* installed size is 76.6Mb sub-directories of 1Mb or more: libs  76.3Mb
    + This package contains compiled stan code, which is much faster, but uses
      these large libs folders.
