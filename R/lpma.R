lambda_range <- function(x, y, vi,
                         standardize = TRUE,
                         intercept = TRUE,
                         lambda.min.ratio = NULL,
                         lambda.max.ratio = 1.1,
                         nlambda = 100,
                         tau2 = NULL) {

  x <- as.matrix(x)
  y <- as.numeric(y)
  vi <- as.numeric(vi)

  nobs <- length(y)
  nvars <- ncol(x)

  if (nrow(x) != nobs || length(vi) != nobs)
    stop("Incompatible dimensions of x, y, and vi.")

  if (any(!is.finite(c(x, y, vi))) || any(vi <= 0))
    stop("Inputs must be finite and vi must be positive.")

  if (nvars < 1L)
    stop("At least one moderator is required.")

  if (nlambda < 2L)
    stop("'nlambda' must be at least 2.")

  ## Default ratio, following glmnet convention
  if (is.null(lambda.min.ratio)) {
    lambda.min.ratio <- if (nobs > nvars) 0.0001 else 0.01
  }

  if (!is.finite(lambda.min.ratio) ||
      lambda.min.ratio <= 0 ||
      lambda.min.ratio >= 1) {
    stop("'lambda.min.ratio' must be between 0 and 1.")
  }

  ## Standardize moderators as in lrma()
  if (standardize) {

    center <- colMeans(x)
    scale <- apply(x, 2, sd)

    if (any(!is.finite(scale)) || any(scale == 0))
      stop("Moderators must have nonzero standard deviations.")

    xs <- sweep(x, 2, center, "-")
    xs <- sweep(xs, 2, scale, "/")

  } else {
    xs <- x
  }

  ## Estimate tau2 under the intercept-only model
  ## using REML, unless supplied.
  if (is.null(tau2)) {

    X0 <- if (intercept) {
      matrix(1, nrow = nobs, ncol = 1)
    } else {
      matrix(numeric(0), nrow = nobs, ncol = 0)
    }

    neg_reml <- function(tau2) {

      w <- 1 / (vi + tau2)

      if (intercept) {
        beta0 <- sum(w * y) / sum(w)
        residuals <- y - beta0
        adjustment <- log(sum(w))
      } else {
        residuals <- y
        adjustment <- 0
      }

      0.5 * (
        sum(log(vi + tau2)) +
          adjustment +
          sum(w * residuals^2)
      )
    }

    ## Include tau2 = 0 as a possible boundary solution
    upper <- max(var(y), mean(vi), 1e-4) * 100

    opt <- optim(
      par = 0,
      fn = function(z) neg_reml(z[1]),
      method = "L-BFGS-B",
      lower = 0,
      upper = upper
    )

    tau2 <- opt$par[1]
  }

  if (length(tau2) != 1L ||
      !is.finite(tau2) || tau2 < 0)
    stop("'tau2' must be a nonnegative finite number.")

  ## Inverse-variance weights
  w <- 1 / (vi + tau2)

  ## Unpenalized intercept
  if (intercept) {
    beta0 <- sum(w * y) / sum(w)
  } else {
    beta0 <- 0
  }

  residuals <- y - beta0

  ## KKT condition for all moderator coefficients = 0
  score <- drop(crossprod(xs, w * residuals))

  lambda.max <- max(abs(score))
  lambda.min <- lambda.max * lambda.min.ratio
  lambda.max <- lambda.max.ratio * lambda.max
  ## Logarithmically spaced penalty sequence
  lambda <- if (lambda.max > 0) {
    exp(seq(
      log(lambda.max),
      log(lambda.min),
      length.out = nlambda
    ))
  } else {
    rep(0, nlambda)
  }

  return(lambda)
}

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
#' `lrma` compute its own lambda sequence.
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
#' @return A `list` object of class `lrma`, with the following structure:
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
#' suppressWarnings({res <- lrma(d~., data = df, iter = 10)})
lrma <- function(x, ...){
  UseMethod("lrma")
}

#' @method lrma formula
#' @export
#' @rdname lrma
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
  formula_intercept <- attr(terms(formula, data = mf), "intercept") == 1L

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
#' @method lrma default
#' @export
#' @rdname lrma
lrma.default <- function(x,
                         y,
                         vi,
                         standardize = TRUE,
                         intercept = TRUE,
                         lambda,
                         tau2_start = NULL,
                         beta_start = NULL,
                         control = list(
                           maxit = 200,
                           tol = 1e-7,
                           cd.maxit = 10000,
                           cd.tol = 1e-8,
                           tau2.upper = NULL
                         ),
                         ...) {

  if (missing(lambda) || length(lambda) != 1L ||
      !is.finite(lambda) || lambda < 0)
    stop("'lambda' must be a finite, nonnegative number.")

  x <- as.matrix(x)
  storage.mode(x) <- "double"
  y <- as.numeric(y)
  vi <- as.numeric(vi)

  n <- length(y)
  p <- ncol(x)

  if (nrow(x) != n || length(vi) != n)
    stop("Incompatible dimensions.")

  if (any(!is.finite(c(x, y, vi))) || any(vi <= 0))
    stop("Invalid data or sampling variances.")

  if (p < 1L)
    stop("At least one moderator is required.")

  if (is.null(colnames(x)))
    colnames(x) <- paste0("x", seq_len(p))

  moderator_names <- colnames(x)

  ## --------------------------------------------------
  ## 1. Standardize moderators
  ## --------------------------------------------------

  center <- rep(0, p)
  scale <- rep(1, p)

  names(center) <- moderator_names
  names(scale) <- moderator_names

  if (standardize) {

    center <- colMeans(x)
    scale <- apply(x, 2, sd)

    if (any(!is.finite(scale)) || any(scale == 0))
      stop("Moderators must have nonzero variance.")

    xs <- sweep(x, 2, center, "-")
    xs <- sweep(xs, 2, scale, "/")

  } else {
    xs <- x
  }

  X <- if (intercept) {
    cbind("(Intercept)" = 1, xs)
  } else {
    xs
  }

  q <- ncol(X)

  ## --------------------------------------------------
  ## 2. Controls
  ## --------------------------------------------------

  defaults <- list(
    maxit = 200,
    tol = 1e-7,
    cd.maxit = 10000,
    cd.tol = 1e-8,
    tau2.upper = NULL
  )

  ctrl <- modifyList(defaults, control)

  ## --------------------------------------------------
  ## 3. Initialize parameters
  ## --------------------------------------------------

  if (is.null(tau2_start)) {
    tau2 <- max(0, var(y) - mean(vi))
  } else {
    tau2 <- tau2_start
  }

  if (!is.finite(tau2) || tau2 < 0)
    stop("Invalid tau2_start.")

  beta <- if (is.null(beta_start)) {
    numeric(q)
  } else {
    as.numeric(beta_start)
  }

  if (length(beta) != q)
    stop("Incorrect length of beta_start.")

  ## --------------------------------------------------
  ## 4. Penalized marginal negative log-likelihood
  ## --------------------------------------------------

  objective <- function(beta, tau2) {

    V <- vi + tau2
    r <- y - drop(X %*% beta)

    0.5 * (
      sum(log(V)) +
        sum(r^2 / V)
    ) +
      lambda * sum(abs(
        if (intercept) beta[-1L] else beta
      ))
  }

  ## --------------------------------------------------
  ## 5. Alternating optimization
  ## --------------------------------------------------

  converged <- FALSE
  history <- numeric(ctrl$maxit)

  for (iter in seq_len(ctrl$maxit)) {

    beta_old <- beta
    tau2_old <- tau2

    ## Step A: coordinate descent for beta
    w <- 1 / (vi + tau2)

    cd <- lasso_cd_cpp(
      X = X,
      y = y,
      w = w,
      lambda = lambda,
      intercept = intercept,
      beta_start = beta,
      tol = ctrl$cd.tol,
      maxit = ctrl$cd.maxit
    )

    beta <- cd$beta

    ## Step B: optimize tau2 conditional on beta
    r <- y - drop(X %*% beta)

    negll_tau2 <- function(t) {
      0.5 * sum(
        log(vi + t) +
          r^2 / (vi + t)
      )
    }

    ## Numerical upper bound
    upper <- ctrl$tau2.upper

    if (is.null(upper)) {
      upper <- max(
        var(y),
        mean(vi),
        tau2,
        1e-4
      ) * 100
    }

    opt <- optim(
      par = tau2,
      fn = function(z) negll_tau2(z[1]),
      method = "L-BFGS-B",
      lower = 0,
      upper = upper
    )

    tau2 <- opt$par[1]

    history[iter] <- objective(beta, tau2)

    ## Convergence check
    beta_change <- max(abs(beta - beta_old))
    tau_change <- abs(tau2 - tau2_old)

    if (max(beta_change, tau_change) <
        ctrl$tol * (1 + max(abs(c(beta_old, tau2_old))))) {

      converged <- TRUE
      break
    }
  }

  history <- history[seq_len(iter)]

  ## --------------------------------------------------
  ## 6. Back-transform coefficients
  ## --------------------------------------------------

  beta_std <- beta
  names(beta_std) <- colnames(X)

  if (standardize) {

    if (intercept) {

      b <- beta[-1L] / scale

      b0 <- beta[1L] - sum(center * b)

      beta_original <- c(
        "(Intercept)" = b0,
        setNames(b, moderator_names)
      )

    } else {

      beta_original <- setNames(
        beta / scale,
        moderator_names
      )
    }

  } else {

    beta_original <- beta_std
  }

  ## --------------------------------------------------
  ## 7. Fitted values and output
  ## --------------------------------------------------

  X_original <- if (intercept) {
    cbind("(Intercept)" = 1, x)
  } else {
    x
  }

  fitted_values <- drop(X_original %*% beta_original)

  out <- list(
    coefficients = beta_original,
    coef_std = beta_std,
    tau2 = tau2,
    lambda = lambda,
    fitted.values = fitted_values,
    residuals = y - fitted_values,
    x = x,
    y = y,
    vi = vi,
    standardize = standardize,
    center = center,
    scale = scale,
    intercept = intercept,
    objective = objective(beta, tau2),
    objective_history = history,
    convergence = if (converged && cd$converged &&
                      opt$convergence == 0L) 0L else 1L,
    iterations = iter,
    cd.converged = cd$converged,
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
                            lambda = NULL,
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

  formula_intercept <- attr(terms(formula, data = mf), "intercept") == 1L

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
  if(is.null(lambda)){
    lambda <- lambda_range(x = X, y = y, vi = vi_vec,
                             standardize = standardize,
                             intercept = intercept,
                             lambda.max.ratio = 1.1,
                             nlambda = 100)
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


## Soft-thresholding operator
soft_threshold <- function(z, lambda) {
  sign(z) * pmax(abs(z) - lambda, 0)
}


## Weighted LASSO coordinate descent
lasso_cd <- function(X, y, w, lambda,
                     intercept = TRUE,
                     beta_start = NULL,
                     tol = 1e-8,
                     maxit = 10000) {

  n <- length(y)
  p <- ncol(X)

  if (is.null(beta_start)) {
    beta <- numeric(p)
  } else {
    beta <- as.numeric(beta_start)
  }

  ## X includes intercept as first column when requested
  penalized <- if (intercept) 2:p else seq_len(p)
  if (intercept && p == 1L)
    penalized <- integer(0)

  ## Maintain residuals for efficient coordinate updates
  r <- y - drop(X %*% beta)

  ## Weighted squared column norms
  z <- colSums(w * X^2)

  converged <- FALSE

  for (iter in seq_len(maxit)) {

    beta_old <- beta

    ## Update unpenalized intercept
    if (intercept) {

      j <- 1L

      r <- r + X[, j] * beta[j]

      beta[j] <- sum(w * r) / z[j]

      r <- r - X[, j] * beta[j]
    }

    ## Update penalized moderator coefficients
    for (j in penalized) {

      r <- r + X[, j] * beta[j]

      rho <- sum(w * X[, j] * r)

      beta[j] <- if (z[j] > 0) {
        soft_threshold(rho, lambda) / z[j]
      } else {
        0
      }

      r <- r - X[, j] * beta[j]
    }

    if (max(abs(beta - beta_old)) <
        tol * (1 + max(abs(beta_old)))) {

      converged <- TRUE
      break
    }
  }

  list(
    beta = beta,
    residuals = r,
    iterations = iter,
    converged = converged
  )
}



predict.lrma <- function(object,
                         newdata = NULL,
                         ...) {
  if (is.null(newdata)) {
    return(object$fitted.values)
  }

  ## Construct moderator matrix
  if (!is.null(object$formula)) {

    ## Formula interface: preserve factor coding,
    ## transformations, and interaction terms.
    tt <- delete.response(terms(object$formula))

    X <- model.matrix(
      tt,
      data = newdata
    )

    if ("(Intercept)" %in% colnames(X)) {
      X <- X[, colnames(X) != "(Intercept)",
             drop = FALSE]
    }

  } else {

    ## Default interface
    X <- as.matrix(newdata)
    storage.mode(X) <- "double"
  }

  ## Match moderators by name
  moderator_names <- if (object$intercept) {
    names(object$coefficients)[-1L]
  } else {
    names(object$coefficients)
  }

  if (!all(moderator_names %in% colnames(X))) {
    stop("newdata is missing required moderators.")
  }

  X <- X[, moderator_names, drop = FALSE]

  ## Construct prediction design matrix
  if (object$intercept) {
    X <- cbind("(Intercept)" = 1, X)
  }

  return(drop(X %*% object$coefficients))
}


predict.cv.lrma <- function(object,
                            newdata = NULL,
                            lambda = "min",
                            ...) {
  ## Resolve lambda
  if (is.character(lambda)) {

    lambda <- match.arg(lambda, c("min", "1se"))

    if (lambda == "min") {
      fit <- object$fit
    } else {
      fit <- object$fit.1se
    }

  } else if (is.numeric(lambda) &&
             length(lambda) == 1L &&
             is.finite(lambda) &&
             lambda >= 0) {

    ## Reuse already fitted models where possible
    if (isTRUE(lambda == object$lambda.min)) {
      fit <- object$fit
    } else if (isTRUE(lambda == object$lambda.1se)) {
      fit <- object$fit.1se
    } else {
      ## Requires training data stored in CV object
      if (is.null(object$fit$x) ||
          is.null(object$fit$y) ||
          is.null(object$fit$vi)) {
        stop(
          "The CV object does not contain training data. ",
          "Store x, y, and vi in cv.lrma.default() ",
          "to enable prediction at arbitrary lambda values."
        )
      }

      fit <- lrma.default(
        x = object$fit$x,
        y = object$fit$y,
        vi = object$fit$vi,
        standardize = object$standardize,
        intercept = object$intercept,
        lambda = lambda,
        ...
      )

      if (fit$convergence != 0L) {
        warning("Refitted model did not converge.")
      }

      ## Preserve formula metadata for prediction
      fit$formula <- object$formula
    }

  } else {

    stop(
      "'lambda' must be 'min', '1se', ",
      "or a nonnegative numeric value."
    )
  }
  return(predict(fit, newdata = newdata))
}
