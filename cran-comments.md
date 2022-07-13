# pema 0.1.2

* Updated maintainer email address
* Use str2lang("pema::function") instead of quote(function) so pema functions
  work without attaching package namespace
* Coerce mean and sd to array so stan accepts cases with 1 moderator
* Throw informative error in cases with 0 moderators

## Test environments
* local x86_64-pc-linux-gnu, R 4.1.2
* GitHub actions, windows-latest (release)
* GitHub actions, macOS-latest (release)
* GitHub actions, ubuntu-20.04 (release)
* GitHub actions, ubuntu-20.04 (devel)
* devtools::check_win_devel()
* devtools::check_win_release()
* devtools::check_win_oldrelease()

## R CMD check results

0 errors | 0 warnings | 1 notes

* New maintainer
  + I updated my email address to reflect my new affiliation
