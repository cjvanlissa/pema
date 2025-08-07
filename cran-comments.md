# pema 0.1.5

* Address Prof. Dr. Brian Ripley's comments / fix CRAN errors caused by
  dependency `metaforest` no longer Import-ing `metafor`

## Test environments

* local Windows 11, R 4.5.0
* devtools::check_win_devel()
* devtools::check_win_release()
* devtools::check_win_oldrelease()
* rhub::rhub_check(platforms = c(1,2,3,4))
* GitHub actions, windows-latest (release)
* GitHub actions, macOS-latest (release)
* GitHub actions, ubuntu-20.04 (release)
* GitHub actions, ubuntu-20.04 (devel)

## R CMD check results

0 errors | 0 warnings | 1 note

* GNU make is a SystemRequirements.
