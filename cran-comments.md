# pema 0.1.4

* Incorporate necessary changes to maintain compatibility with rstan
* Make Suggested packages truly optional, add test without Suggests
* Add tutorial vignette

## Test environments
* local Windows 11, R 4.4.3
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
