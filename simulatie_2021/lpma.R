#' Conduct LASSO Penalized Meta-Analysis
#'
#' This function conducts LASSO penalized meta-regression (Van Lissa, in prep).
#' L1-penalization is used to shrink the regression coefficients of irrelevant
#' moderators towards zero, within the context of restricted maximum likelihood
#' estimation.
#' @param formula An object of class `formula` (or one that can be coerced to
#' that class), see \code{\link[stats]{lm}}.
#' @param data Either a `data.frame` containing the variables in the model,
#' see \code{\link[stats]{lm}}, or a `list` of multiple imputed `data.frame`s,
#' or an object returned by \code{\link[mice]{mice}}.
#' @param vi Character. Name of the column in the \code{data} that
#' contains the variances of the effect sizes. This column will be removed from
#' the data prior to analysis. Defaults to \code{"vi"}.
# @param study Character. Name of the column in the
# \code{data} that contains the study id. Use this when the data includes
# multiple effect sizes per study. This column can be a vector of integers, or
# a factor. This column will be removed from the data prior to analysis.
# See \code{Details} for more information about analyzing dependent data.
#' @param standardize Either a logical argument or a list. If `standardize` is
#' logical, it controls whether all predictors are standardized prior to
#' analysis or not. Parameter estimates are restored to the predictors' original
#' scale. Alternatively, users can provide a list to `standardize` to gain
#' more control over the standardization process. In this case, it is assumed
#' that the standardization has already taken place. This list must have two
#' elements: `list(center = c(mean(X1)
#' , mean(X2), mean(X...)), scale = c(sd(X1), sd(X2), sd(X...)))`. It is used
#' only to restore parameter estimates to the original scale of the predictors.
#' This is useful, e.g., to standardize continuous and dichotomous variables
#' separately. Defaults to \code{TRUE}, which is recommended so that shrinking
#' affects all parameters similarly.
#' @param intercept Logical, indicating whether or not an intercept should be included
#' in the model.
#' @param lambda A user supplied lambda sequence. Typical usage is to have
#' `lpma` compute its own lambda sequence.
#' @param ... Additional arguments passed on to other functions.
# @references
# Van Lissa, C. J., van Erp, S., & Clapper, E. B. (2023). Selecting relevant
# moderators with Bayesian regularized meta-regression. Research Synthesis
# Methods. \doi{10.31234/osf.io/6phs5}
#
# Park, T., & Casella, G. (2008). The Bayesian Lasso. Journal of the American
# Statistical Association, 103(482), 681–686. \doi{10.1198/016214508000000337}
#
# Carvalho, C. M., Polson, N. G., & Scott, J. G. (2010). The horseshoe
# estimator for sparse signals. Biometrika, 97(2), 465–480.
# \doi{10.1093/biomet/asq017}
#
# Piironen, J., & Vehtari, A. (2017). Sparsity information and regularization
# in the horseshoe and other shrinkage priors. Electronic Journal of
# Statistics, 11(2). \url{https://projecteuclid.org/journals/electronic-journal-of-statistics/volume-11/issue-2/Sparsity-nformation-and-regularization-in-the-horseshoe-and-other-shrinkage/10.1214/17-EJS1337SI.pdf}
#' @return A `list` object of class `lpma`, with the following structure:
#'   out <- list(
# tau2 = tau2_hat,
# tau = sqrt(tau2_hat),
# lambda = lambda,
# fitted.values = fitted,
# residuals = residuals,
# x = x,
# y = y,
# vi = vi,
# standardize = standardize,
# center = center,
# scale = scale,
# intercept = intercept,
# objective = opt$value,
# convergence = opt$convergence,
# message = opt$message,
# optimizer = opt,
# call = match.call()
#' ```
#' list(
#'   coefficients # A numeric matrix with parameter estimates; these are
#'                # interpreted as regression coefficients, except tau2 and tau,
#'                # which are interpreted as the residual variance and standard
#'                # deviation, respectively.
#'   coef_std     # A numeric matrix with standardized parameter estimates;
#'                # these are interpreted as standardized regression
#'                # coefficients, except tau2 and tau, which are interpreted as
#'                # the residual variance and standard deviation, respectively.
#'   tau2         # The estimated residual variance.
#'   formula      # The formula used to estimate the model
#'   terms        # The predictor terms in the formula
#'   X            # Numeric matrix of moderator variables
#'   Y            # Numeric vector with effect sizes
#'   vi           # Numeric vector with effect size variances
#'   tau2         # Numeric, estimated tau2
#'   R2           # Numeric, estimated heterogeneity explained by the moderators
#'   k            # Numeric, number of effect sizes
#'   study        # Numeric vector with study id numbers
#' )
#' ```
#' @export
#' @examples
#' data("curry")
#' df <- curry[c(1:5, 50:55), c("d", "vi", "sex", "age", "donorcode")]
#' suppressWarnings({res <- lpma(d~., data = df, iter = 10)})
lpma <- function(x, ...){
  UseMethod("lpma")
}

#' @method lpma formula
#' @export
#' @rdname lpma
lrma.formula <- function(formula,
                         data,
                         vi = "vi",
                         standardize = TRUE,
                         intercept = TRUE,
                         lambda,
                         ...) {

  if (missing(data)) {
    stop("'data' must be supplied.")
  }

  ## Construct model frame for outcome + moderators
  mf <- model.frame(
    formula = formula,
    data = data,
    na.action = na.pass
  )

  y <- model.response(mf)

  ## Construct design matrix from formula
  X <- model.matrix(formula, data = mf)

  ## Determine whether formula contains an intercept
  formula_intercept <- attr(terms(formula), "intercept") == 1L

  ## Remove model.matrix intercept because the default method
  ## handles the intercept separately.
  if (formula_intercept && "(Intercept)" %in% colnames(X)) {
    X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
  }

  ## If formula explicitly says ~ 0 + ..., respect that unless
  ## intercept was supplied explicitly.
  if (!formula_intercept) {
    intercept <- FALSE
  }

  ## Sampling variances can be specified by column name
  if (is.character(vi) && length(vi) == 1L) {
    if (!vi %in% names(data)) {
      stop("Column '", vi, "' not found in 'data'.")
    }
    vi_vec <- data[[vi]]
  } else {
    vi_vec <- vi
  }

  ## Complete-case filtering across yi, vi, and moderators
  keep <- complete.cases(y, vi_vec, X)

  if (!all(keep)) {
    warning(sum(!keep), " observation(s) removed due to missing values.")
  }

  fit <- lrma.default(
    x = X[keep, , drop = FALSE],
    y = y[keep],
    vi = vi_vec[keep],
    standardize = standardize,
    intercept = intercept,
    lambda = lambda,
    ...
  )

  fit$call <- match.call()
  fit$formula <- formula
  fit$na.action <- which(!keep)
  return(fit)
}


#' @param x An k x m numeric matrix, where k is the number of effect sizes and m
#' is the number of moderators.
#' @param y A numeric vector of k effect sizes.
#' @param tau2_start Starting value for `tau2`. Defaults to the effect size
#' variance minus the mean sampling variance, with a minimum of `1e-8`.
#' @param beta_start Vector of starting values for `coefficients`. Defaults to
#' the random effects-weighted regression coefficients based on `tau2_start`.
#' @method lpma default
#' @export
#' @rdname lpma
lrma.default <- function(x,
                         y,
                         vi,
                         standardize = TRUE,
                         intercept = TRUE,
                         lambda,
                         tau2_start = NULL,
                         beta_start = NULL,
                         control = list(maxit = 10000),
                         ...) {

  if (missing(lambda)) {
    stop("'lambda' must be supplied. Use cv.lrma() to tune lambda.")
  }

  if (length(lambda) != 1L ||
      !is.finite(lambda) ||
      lambda < 0) {
    stop("'lambda' must be a single non-negative number.")
  }

  y <- as.numeric(y)
  vi <- as.numeric(vi)

  if (is.data.frame(x)) {
    x <- data.matrix(x)
  } else {
    x <- as.matrix(x)
  }

  n <- length(y)

  if (length(vi) != n) {
    stop("'y' and 'vi' must have the same length.")
  }

  if (nrow(x) != n) {
    stop("'x' must have one row per effect size.")
  }

  if (any(!is.finite(y))) {
    stop("'y' must contain only finite values.")
  }

  if (any(!is.finite(vi)) || any(vi <= 0)) {
    stop("'vi' must contain finite, positive sampling variances.")
  }

  if (any(!is.finite(x))) {
    stop("'x' must contain only finite values.")
  }

  if (ncol(x) == 0L) {
    stop("At least one moderator must be supplied.")
  }

  if (is.null(colnames(x))) {
    colnames(x) <- paste0("x", seq_len(ncol(x)))
  }

  moderator_names <- colnames(x)

  # Standardize moderators
  center <- rep(0, ncol(x))
  scale <- rep(1, ncol(x))

  names(center) <- moderator_names
  names(scale) <- moderator_names

  if (standardize) {

    center <- colMeans(x)
    scale <- apply(x, 2L, sd)

    if (any(!is.finite(scale)) || any(scale == 0)) {
      stop("All moderators must have non-zero variance when ",
           "'standardize = TRUE'.")
    }

    x_fit <- sweep(x, 2L, center, "-")
    x_fit <- sweep(x_fit, 2L, scale, "/")

  } else {

    x_fit <- x
  }

  # Design matrix
  if (intercept) {
    X <- cbind("(Intercept)" = 1, x_fit)
  } else {
    X <- x_fit
  }

  p <- ncol(X)

  if (n <= p) {
    warning(
      "The number of observations is not larger than the number of fixed-effect coefficients. REML estimation may be unstable."
    )
  }

  # tau² starting value
  if (is.null(tau2_start)) {
    tau2_start <- max(
      1e-8,
      var(y) - mean(vi)
    )
  }

  if (!is.finite(tau2_start) || tau2_start <= 0) {
    tau2_start <- 1e-4
  }

  # beta starting values
  if (is.null(beta_start)) {

    w0 <- 1 / (vi + tau2_start)

    XtWX <- crossprod(X, X * w0)
    XtWy <- crossprod(X, y * w0)

    beta_start <- tryCatch(
      drop(solve(XtWX, XtWy)),
      error = function(e) rep(0, p)
    )
  }

  if (length(beta_start) != p) {
    stop("'beta_start' must have length ", p, ".")
  }

  # Penalized REML objective function
  objective <- function(par) {

    beta <- par[seq_len(p)]
    tau2 <- exp(par[p + 1L])

    Vdiag <- vi + tau2
    w <- 1 / Vdiag

    residuals <- y - drop(X %*% beta)

    # Gaussian likelihood component
    obj <- 0.5 * (
      sum(log(Vdiag)) +
        sum(w * residuals^2)
    )

    # REML adjustment
    XtWX <- crossprod(X, X * w)

    R <- tryCatch(
      chol(XtWX),
      error = function(e) NULL
    )

    if (is.null(R)) {
      return(.Machine$double.xmax / 100)
    }

    logdet_XtWX <- 2 * sum(log(diag(R)))

    obj <- obj + 0.5 * logdet_XtWX

    # LASSO penalty (ignoring the intercept)
    if (intercept) {
      penalty_beta <- beta[-1L]
    } else {
      penalty_beta <- beta
    }

    obj + lambda * sum(abs(penalty_beta))
  }

  par_start <- c(beta_start, log(tau2_start))
  # Estimate model
  opt <- optim(
    par = par_start,
    fn = objective,
    method = "Nelder-Mead",
    control = control,
    ...
  )

  beta_std <- opt$par[seq_len(p)]
  tau2_hat <- exp(opt$par[p + 1L])

  names(beta_std) <- colnames(X)

  # Restore coefficients to original scale
  if (standardize) {
    if (intercept) {
      beta_mod_std <- beta_std[-1L]
      beta_mod <- beta_mod_std / scale
      beta_intercept <- beta_std[1L] - sum(beta_mod_std * center / scale)
      beta <- c(
        "(Intercept)" = beta_intercept,
        beta_mod
      )
    } else {
      beta <- beta_std / scale
    }
    names(beta)[if (intercept) -1L else seq_along(beta)] <- moderator_names
  } else {
    beta <- beta_std
  }

  # Predictions on original scale
  if (intercept) {
    X_original <- cbind("(Intercept)" = 1, x)
  } else {
    X_original <- x
  }
  fitted <- drop(X_original %*% beta)
  residuals <- y - fitted

  out <- list(
    coefficients = beta,
    coef_std = beta_std,
    tau2 = tau2_hat,
    tau = sqrt(tau2_hat),
    lambda = lambda,
    fitted.values = fitted,
    residuals = residuals,
    x = x,
    y = y,
    vi = vi,
    standardize = standardize,
    center = center,
    scale = scale,
    intercept = intercept,
    objective = opt$value,
    convergence = opt$convergence,
    message = opt$message,
    optimizer = opt,
    call = match.call()
  )

  class(out) <- c("lrma", class(out))
  return(out)
}



coef.lrma <- function(object, standardized = FALSE, ...) {

  if (standardized) {
    object$coef_std
  } else {
    object$coefficients
  }
}


fitted.lrma <- function(object, ...) {
  object$fitted.values
}


residuals.lrma <- function(object, ...) {
  object$residuals
}


print.lrma <- function(x, ...) {

  cat("LASSO-regularized meta-regression\n\n")

  cat("lambda:", signif(x$lambda, 5), "\n")
  cat("tau^2:  ", signif(x$tau2, 5), "\n")
  cat("tau:    ", signif(x$tau, 5), "\n\n")

  cat("Coefficients:\n")
  print(x$coefficients)

  if (x$convergence != 0L) {
    cat("\nOptimizer convergence code:", x$convergence, "\n")
  }

  invisible(x)
}


## ============================================================
## K-fold cross-validation for lrma()
## ============================================================

cv.lrma <- function(x, ...) {
  UseMethod("cv.lrma")
}


## ------------------------------------------------------------
## Formula interface
## ------------------------------------------------------------

cv.lrma.formula <- function(formula,
                            data,
                            vi = "vi",
                            standardize = TRUE,
                            intercept = TRUE,
                            lambda = exp(seq(log(1e-4),
                                             log(10),
                                             length.out = 50)),
                            folds = 5,
                            seed = NULL,
                            metric = c("weighted_mse", "mse"),
                            ...) {

  metric <- match.arg(metric)

  mf <- model.frame(
    formula = formula,
    data = data,
    na.action = na.pass
  )

  y <- model.response(mf)
  X <- model.matrix(formula, data = mf)

  formula_intercept <- attr(terms(formula), "intercept") == 1L

  if (formula_intercept && "(Intercept)" %in% colnames(X)) {
    X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
  }

  if (!formula_intercept) {
    intercept <- FALSE
  }

  if (is.character(vi) && length(vi) == 1L) {

    if (!vi %in% names(data)) {
      stop("Variable '", vi, "' not found in 'data'.")
    }

    vi_vec <- data[[vi]]

  } else {

    vi_vec <- vi
  }

  keep <- complete.cases(y, vi_vec, X)

  if (!all(keep)) {
    warning(sum(!keep), " observation(s) removed because of missing values.")
  }

  out <- cv.lrma.default(
    x = X[keep, , drop = FALSE],
    y = y[keep],
    vi = vi_vec[keep],
    standardize = standardize,
    intercept = intercept,
    lambda = lambda,
    folds = folds,
    seed = seed,
    metric = metric,
    ...
  )

  out$call <- match.call()
  out$formula <- formula
  out$na.action <- which(!keep)

  out
}


## ------------------------------------------------------------
## Default interface
## ------------------------------------------------------------

cv.lrma.default <- function(x,
                            y,
                            vi,
                            standardize = TRUE,
                            intercept = TRUE,
                            lambda = exp(seq(log(1e-4),
                                             log(10),
                                             length.out = 50)),
                            folds = 5,
                            seed = NULL,
                            metric = c("weighted_mse", "mse"),
                            ...) {

  metric <- match.arg(metric)

  y <- as.numeric(y)
  vi <- as.numeric(vi)

  if (is.data.frame(x)) {
    x <- data.matrix(x)
  } else {
    x <- as.matrix(x)
  }

  n <- length(y)

  if (length(vi) != n) {
    stop("'y' and 'vi' must have the same length.")
  }

  if (nrow(x) != n) {
    stop("'x' must have one row per effect size.")
  }

  if (folds < 2L || folds > n) {
    stop("'folds' must be between 2 and the number of observations.")
  }

  lambda <- sort(unique(as.numeric(lambda)))

  if (any(!is.finite(lambda)) || any(lambda < 0)) {
    stop("All lambda values must be finite and non-negative.")
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }

  ## Balanced random folds
  fold_id <- sample(
    rep(seq_len(folds), length.out = n)
  )

  fold_error <- matrix(
    NA_real_,
    nrow = length(lambda),
    ncol = folds
  )

  ## ----------------------------------------------------------
  ## Cross-validation
  ##
  ## Importantly, standardization happens INSIDE lrma() using
  ## the training data only. Thus no validation-fold information
  ## leaks into the standardization step.
  ## ----------------------------------------------------------

  for (j in seq_along(lambda)) {

    for (f in seq_len(folds)) {

      test <- fold_id == f
      train <- !test

      fit <- tryCatch(
        lrma.default(
          x = x[train, , drop = FALSE],
          y = y[train],
          vi = vi[train],
          standardize = standardize,
          intercept = intercept,
          lambda = lambda[j],
          ...
        ),
        error = function(e) NULL
      )

      if (is.null(fit) || fit$convergence != 0L) {
        next
      }

      ## Coefficients returned by lrma() are already transformed
      ## back to the original moderator scale.
      if (intercept) {

        X_test <- cbind(
          "(Intercept)" = 1,
          x[test, , drop = FALSE]
        )

      } else {

        X_test <- x[test, , drop = FALSE]
      }

      pred <- drop(
        X_test %*% fit$coefficients
      )

      sq_error <- (y[test] - pred)^2

      if (metric == "mse") {

        fold_error[j, f] <- mean(sq_error)

      } else {

        w <- 1 / vi[test]

        fold_error[j, f] <-
          sum(w * sq_error) / sum(w)
      }
    }
  }

  mean_error <- apply(
    fold_error,
    1L,
    function(z) {
      if (all(!is.finite(z))) {
        Inf
      } else {
        mean(z[is.finite(z)])
      }
    }
  )

  se_error <- apply(
    fold_error,
    1L,
    function(z) {

      z <- z[is.finite(z)]

      if (length(z) <= 1L) {
        return(NA_real_)
      }

      sd(z) / sqrt(length(z))
    }
  )

  min_id <- which.min(mean_error)
  lambda_min <- lambda[min_id]

  ## 1-SE rule
  cutoff <- mean_error[min_id] + se_error[min_id]

  if (is.finite(cutoff)) {

    eligible <- which(mean_error <= cutoff)
    lambda_1se <- max(lambda[eligible])

  } else {

    lambda_1se <- lambda_min
  }

  ## Fit selected models to complete data
  fit_min <- lrma.default(
    x = x,
    y = y,
    vi = vi,
    standardize = standardize,
    intercept = intercept,
    lambda = lambda_min,
    ...
  )

  fit_1se <- lrma.default(
    x = x,
    y = y,
    vi = vi,
    standardize = standardize,
    intercept = intercept,
    lambda = lambda_1se,
    ...
  )

  cv_table <- data.frame(
    lambda = lambda,
    error = mean_error,
    se = se_error
  )

  out <- list(
    lambda = lambda,
    lambda.min = lambda_min,
    lambda.1se = lambda_1se,
    cv = cv_table,
    fold.error = fold_error,
    fold.id = fold_id,
    folds = folds,
    metric = metric,
    standardize = standardize,
    intercept = intercept,
    fit = fit_min,
    fit.1se = fit_1se,
    call = match.call()
  )

  class(out) <- "cv.lrma"

  out
}


## ------------------------------------------------------------
## Print method for CV object
## ------------------------------------------------------------

print.cv.lrma <- function(x, ...) {

  cat("K-fold cross-validation for LASSO meta-regression\n\n")

  cat("Number of folds:", x$folds, "\n")
  cat("CV criterion:    ", x$metric, "\n")
  cat("lambda.min:      ", signif(x$lambda.min, 5), "\n")
  cat("lambda.1se:      ", signif(x$lambda.1se, 5), "\n")

  invisible(x)
}


## ------------------------------------------------------------
## Plot method
## ------------------------------------------------------------

plot.cv.lrma <- function(x, ...) {

  plot(
    log(x$cv$lambda),
    x$cv$error,
    type = "b",
    xlab = "log(lambda)",
    ylab = "Cross-validated prediction error",
    ...
  )

  abline(
    v = log(x$lambda.min),
    lty = 2
  )

  abline(
    v = log(x$lambda.1se),
    lty = 3
  )

  invisible(x)
}
# library(metafor)
# ### calculate log risk ratios and corresponding sampling variances
# dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)
#
# ### fit a mixed-effects model with two moderators (absolute latitude and publication year)
# res_metafor <- rma(yi, vi, mods=cbind(ablat, year), data=dat)
#
# res_rma <- reml_meta_reg(dat$yi, dat$vi, mods = dat[ , c("ablat", "year")], method = "BFGS")
# res_lasso <- lasso_reml_meta_reg(dat$yi, dat$vi, mods = dat[ , c("ablat", "year")], method = "BFGS")
# res_lasso <- cv_lasso_reml_meta_reg(dat$yi, dat$vi, mods = dat[ , c("ablat", "year")])
# res_lasso$fit$coefficients
# res_rma$coefficients$estimate
#
# res_lrma <- cv.lrma(yi ~ ablat + year, data = dat)
# res_brma <- pema::brma(yi ~ ablat + year, data = dat)
