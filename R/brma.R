#' Conduct Bayesian Regularized Meta-Analysis
#'
#' This function uses Bayesian estimation via the \code{stan} function
#' \code{\link[stan]{sampling}} to fit a meta-analytic mixed-effects model with
#' moderators. A lasso or horseshoe prior is used to shrink the regression
#' coefficients of irrrelevant moderators towards zero. See Details.
#' @param formula An object of class `formula` (or one that can be coerced to
#' that class), see \code{\link[stats]{lm}}.
#' @param data Optional data.frame containing the variables in the model, see
#' \code{\link[stats]{lm}}.
#' @param vi Character. Name of the column in the \code{data} that
#' contains the variances of the effect sizes. This column will be removed from
#' the data prior to analysis. Defaults to \code{"vi"}.
#' @param study Character. Name of the column in the
#' \code{data} that contains the study id. Use this when the data includes
#' multiple effect sizes per study. This column can be a vector of integers, or
#' a factor. This column will be removed from the data prior to analysis.
#' See \code{Details} for more information about analyzing dependent data.
#' @param method Character, indicating the type of regularizing prior to use.
#' Supports one of \code{c("lasso", "hs")}, see Details. Defaults to
#' \code{"lasso"}.
#' @param scale Logical, indicating whether or not to scale the predictors
#' (defaults to \code{TRUE}, which is recommended so that shrinking affects all
#' parameters similarly.
#' @param prior Numeric vector, specifying the prior to use. Note that the
#' different \code{method}s require this vector to contain specific named
#' elements.
# @param iter A positive integer specifying the number of iterations for each
# chain (including warmup). Defaults to 2000.
#  the model statement. Defaults to .5.
# @param chains A positive integer specifying the number of Markov chains.
# Defaults to 4.
#' @param ... Additional arguments passed on to [rstan::sampling()].
#' Use this, e.g., to override default arguments of that function.
#' @details The following models are available:
#' \describe{
#'   \item{lasso}{ More info about the lasso model}
#'   \item{hs}{ More info about the horseshoe model}
#' }
#' @export
#' @examples
#' set.seed(8)
#' SimulateSMD()
#' SimulateSMD(k_train = 50, distribution = "bernoulli")
#' SimulateSMD(distribution = "bernoulli", model = es * x[ ,1] * x[ ,2])
brma <-
  function(formula,
           data,
           vi = "vi",
           study = "study",
           method = "hs",
           scale = TRUE,
           prior = switch(method, "lasso" = c(lasso_df = 1, lasso_scale = 1),
                                  "hs" = c(hs_df = 1, hs_df_global = 1, hs_df_slab = 4, hs_scale_global = 1, hs_scale_slab = 1)),
           ...) {
    # real<lower=0> hs_df;  // local degrees of freedom
    # real<lower=0> hs_df_global;  // global degrees of freedom
    # real<lower=0> hs_df_slab;  // slab degrees of freedom
    # real<lower=0> hs_scale_global;  // global prior scale
    # real<lower=0> hs_scale_slab;  // slab prior scale
  mf <- match.call(expand.dots = FALSE)
  mf <- mf[c(1L, match(c("formula", "data", "subset", "na.action"), names(mf), 0L))]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  Y <- mf[[1]]
  X <- mf[,-1, drop = FALSE]
  # Add intercept
  #X <- cbind(1, X)
  if(inherits(vi, "character")){
    X[[vi]] <- NULL
    vi <- mf[[vi]]
  }
  if(inherits(study, "character")){
    X[[study]] <- NULL
    study <- mf[[study]]
  }
  se <- sqrt(vi)
  N <- length(Y)
  if(isTRUE(scale)){
    X <- scale(X)
  }

  standat <- c(
    list(
      N = N,
      Y = Y,
      se = se,
      K = ncol(X),
      X = X),
    as.list(prior),
    list(
      N_1 = 20,
      M_1 = 1,
      J_1 = 1:N,
      Z_1_1 = rep(1, N),
      prior_only = FALSE
    )
  )
  cl <- do.call("call",
                c(list(name = "sampling",
                     object = stanmodels[[c(lasso = "lasso_MA", hs = "horseshoe_MA")[method]]],
                     data = standat
                     ),
                  list(...)))
  fit <- eval(cl)
  attr(fit, "type") <- "brma"
  return(fit)
}

