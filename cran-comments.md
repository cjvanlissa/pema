# pema 0.1.4

* Make Suggested packages truly optional, add test without Suggests
* Add vignette
* Maintain compatibility with rstan

## Test environments
* local Windows 11, R 4.2.2
* local x86_64-pc-linux-gnu, R 4.1.2
* GitHub actions, windows-latest (release)
* GitHub actions, macOS-latest (release)
* GitHub actions, ubuntu-20.04 (release)
* GitHub actions, ubuntu-20.04 (devel)
* devtools::check_win_devel()
* devtools::check_win_release()
* devtools::check_win_oldrelease()
* rhub::check_for_cran()

## R CMD check results

0 errors | 0 warnings | 2 notes

* GNU make is a SystemRequirements.
* Some C++ header files contain CR or CRLF line endings.
