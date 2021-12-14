#' pema: Conduct penalized meta-regression.
#'
#' @description Penalized meta-regression shrinks the regression slopes of
#' irrelevant moderators towards zero (Van Lissa & Van Erp, 2021).
#'
#' @docType package
#' @name pema-package
#' @aliases pema
#' @useDynLib pema, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @importFrom rstan sampling
#'
#' @references
#' Van Lissa, C. J., & van Erp, S. (2021, December 9). Select relevant
#' moderators using Bayesian regularized meta-regression.
#' \doi{10.31234/osf.io/6phs5}
#'
#' Stan Development Team (NA). RStan: the R interface to Stan. R package version
#' 2.26.2. \url{https://mc-stan.org}
#'
NULL
