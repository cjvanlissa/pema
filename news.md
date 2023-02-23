# pema 0.1.3

* Updated reference to published validation paper

# pema 0.1.2

* Updated maintainer email address
* Use str2lang("pema::function") instead of quote(function) so pema functions
  work without attaching package namespace
* Coerce mean and sd to array so stan accepts cases with 1 moderator
* Throw informative error in cases with 0 moderators

# pema 0.1.1

* Refactor STAN code
* Support optional intercept
* Support three-level meta-analysis
* Add vignette
* Add Shiny app for visualizing priors
* Add support for multiple imputation

# pema 0.1.0

* First submission to CRAN
