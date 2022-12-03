#' Conduct LASSO penalized Meta-Analysis
#'
#' This function conducts Bayesian regularized meta-regression (Van Lissa & Van
#' Erp, 2021). It uses the \code{stan} function
#' [rstan::sampling] to fit the model. A lasso or horseshoe prior is used to
#' shrink the regression coefficients of irrelevant moderators towards zero.
#' See Details.
#' @param formula An object of class `formula` (or one that can be coerced to
#' that class), see \code{\link[stats]{lm}}.
#' @param data Either a `data.frame` containing the variables in the model,
#' see \code{\link[stats]{lm}}, or a `list` of multiple imputed `data.frame`s,
#' or an object returned by \code{\link[mice]{mice}}.
#' @param vi Character. Name of the column in the \code{data} that
#' contains the variances of the effect sizes. This column will be removed from
#' the data prior to analysis. Defaults to \code{"vi"}.
#' @param study Character. Name of the column in the
#' \code{data} that contains the study id. Use this when the data includes
#' multiple effect sizes per study. This column can be a vector of integers, or
#' a factor. This column will be removed from the data prior to analysis.
#' See \code{Details} for more information about analyzing dependent data.
#' @param ntau Numeric, the numer of `tau2` values to try during
#' cross-validation.
# @details The Bayesian regularized meta-analysis algorithm (Van Lissa & Van
# Erp, 2021) penalizes meta-regression coefficients either via the
# lasso prior (Park & Casella, 2008) or the regularized horseshoe prior
# (Piironen & Vehtari, 2017).
# \describe{
#   \item{lasso}{ The Bayesian equivalent of the lasso penalty is obtained when
#   placing independent Laplace (i.e., double exponential) priors on the
#   regression coefficients centered around zero. The scale of the Laplace
#   priors is determined by a global scale parameter \code{scale}, which
#   defaults to 1 and an inverse-tuning parameter \eqn{\frac{1}{\lambda}}
#   which is given a chi-square prior governed by a degrees of freedom
#   parameter \code{df} (defaults to 1). If \code{standardize = TRUE},
#   shrinkage will
#   affect all coefficients equally and it is not necessary to adapt the
#   \code{scale} parameter. Increasing the \code{df} parameter will allow
#   larger values for the inverse-tuning parameter, leading to less shrinkage.}
#   \item{hs}{ One issue with the lasso prior is that it has relatively light
#   tails. As a result, not only does the lasso have the desirable behavior of
#   pulling small coefficients to zero, it also results in too much shrinkage
#   of large coefficients. An alternative prior that improves upon this
#   shrinkage pattern is the horseshoe prior (Carvalho, Polson & Scott, 2010).
#   The horseshoe prior has an infinitely large spike at zero, thereby pulling
#   small coefficients toward zero but in addition has fat tails, which allow
#   substantial coefficients to escape the shrinkage. The regularized horseshoe
#   is an extension of the horseshoe prior that allows the inclusion of prior
#   information regarding the number of relevant predictors and can
#   be more numerically stable in certain cases (Piironen & Vehtari, 2017).
#   The regularized horseshoe has a global shrinkage parameter that influences
#   all coefficients similarly and local shrinkage parameters that enable
#   flexible shrinkage patterns for each coefficient separately. The local
#   shrinkage parameters are given a Student's t prior with a default \code{df}
#   parameter of 1. Larger values for \code{df} result in lighter tails and
#   a prior that is no longer strictly a horseshoe prior. However, increasing
#   \code{df} slightly might be necessary to avoid divergent transitions in
#   Stan (see also \url{https://mc-stan.org/misc/warnings.html}). Similarly,
#   the degrees of freedom for the Student's t prior on the global shrinkage
#   parameter \code{df_global} can be increased from the default of 1 to, for
#   example, 3 if divergent transitions occur although the resulting
#   prior is then strictly no longer a horseshoe. The scale for the Student's t
#   prior on the global shrinkage parameter \code{scale_global} defaults to 1
#   and can be decreased to achieve more shrinkage. Moreover, if prior
#   information regarding the number of relevant moderators is available, it is
#   recommended to include this information via the \code{relevant_pars}
#   argument by setting it to the expected number of relevant moderators. When
#   \code{relevant_pars} is specified, \code{scale_global} is ignored and
#   instead based on the available prior information. Contrary to the horseshoe
#   prior, the regularized horseshoe applies additional regularization on large
#   coefficients which is governed by a Student's t prior with a
#   \code{scale_slab} defaulting to 2 and \code{df_slab} defaulting to 4.
#   This additional regularization ensures at least some shrinkage of large
#   coefficients to avoid any sampling problems.}
# }
#' @return A `list` object of class `pma`, with the following structure:
#' ```
#' list(
#'   rma          # A random effects meta-analysis of class rma
#'   lasso        # A LASSO-penalized regression model of class cv.glmnet
#'   tau2_cv      # The residual heterogeneity tau2 based on cross-validation
#' )
#' ```
#' @export
# @examples
# data("curry")
# df <- curry[c(1:5, 50:55), c("d", "vi", "sex", "age", "donorcode")]
# suppressWarnings({res <- pma(d~., data = df, ntau = 2)})
#' @importMethodsFrom rstan summary
#' @importFrom stats model.matrix na.omit quantile sd
#' @importFrom metafor rma
#' @importFrom glmnet cv.glmnet
# The line above is just to avoid CRAN warnings that RcppParallel is not
# imported from, despite RcppParallel being a necessary dependency of rstan.
pma <- function(x, ...){
  UseMethod("pma")
}

#' @method pma formula
#' @export
#' @rdname pma
pma.formula <-
  function(formula,
           data,
           vi = "vi",
           study = NULL,
           ntau = 10,
           ...) {
    cl <- match.call()
    # Check for complete data
    if(anyNA(data)) stop("The function pma() requires complete data.")
    # Bookkeeping for columns that should not be in X or Y
    data_cleaned <- .pema_prep_dat(data, vi, study)
    # Make model matrix
    mf <- match.call(expand.dots = FALSE)
    mf <- mf[c(1L, match(c("formula", "subset", "na.action"), names(mf), nomatch = 0L))]
    mf[["data"]] <- data_cleaned$data
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- str2lang("stats::model.frame")
    mf <- eval(mf, parent.frame())
    Y <- mf[[1]]
    mt <- attr(mf, "terms")
    X <- model.matrix(mt, mf)
    if(all(X[,1] == 1)){
      intercept <- TRUE
      X <- X[, -1, drop = FALSE]
    } else {
      intercept <- FALSE
    }
    cl[names(cl) %in% c("formula", "data")] <- NULL
    cl[[1L]] <- str2lang("pema::pma")
    cl[["x"]] <- X
    cl[["y"]] <- Y
    cl[c("vi", "study")] <- data_cleaned[c("vi", "study")]
    cl[["intercept"]] <- intercept
    cl[["formula"]] <- formula
    eval.parent(cl)
}

#' @param x An k x m numeric matrix, where k is the number of effect sizes and m
#' is the number of moderators.
#' @param y A numeric vector of k effect sizes.
#' @param intercept Logical, indicating whether or not an intercept should be included
#' in the model.
#' @method pma default
#' @export
#' @rdname pma
pma.default <-
  function(x,
           y,
           vi,
           study = NULL,
           standardize,
           intercept,
           ntau = 10,
           ...) {
    X <- x
    Y <- y
    # Check for complete data
    if(anyNA(X) | anyNA(Y)) stop("The function pma() requires complete data.")
    # Determine CV folds
    if(!hasArg(foldid)){
      if(!hasArg(nfolds)){
        nfolds <- 10
      }
      if(is.null(study)){
        foldid <- sample(rep(seq(nfolds), length = length(y)))
      } else {
        nstudies <- length(unique(study))
        if(nstudies < nfolds*20) message("Number of unique studies is small relative to the number of cross-validation folds. Decrease 'nfolds'.")
        foldid <- sample(rep(seq(nfolds), length = nstudies))
        foldid <- foldid[study]
      }
    }
    mf_nomod <- metafor::rma(y, vi)
    tau2_max <- mf_nomod$tau2
    mf_allmod <- metafor::rma(yi = y, vi = vi, mods = x)
    tau2_min <- mf_allmod$tau2
    tau_seq <- seq(from = tau2_min, to = tau2_max, length.out = ntau)
    wts <- 1/(vi + tau_seq[1])
    res <- vector(mode = "list", length = ntau)
    res[[1]] <- glmnet::cv.glmnet(x = X, y = y, weights = wts, keep = TRUE, foldid = foldid)
    for(i in 2:ntau){
      wts <- 1/(vi + tau_seq[i])
      res[[i]] <- cv.glmnet(x = X, y = y, weights = wts, foldid = foldid)
    }
    err_min <- sapply(res, function(i){min(i$cvm)})
    tau2_cv <- tau_seq[which.min(err_min)]
    out <- res[[which.min(err_min)]]
    browser()
    out <- c(out,
             list(
               rma = mf_nomod,
               tau2_cv = tau2_cv,
               x = x,
               y = y,
               vi = vi,
               study = study,
               foldid = foldid
               ))
    class(out) <- c("pma", "cv.glmnet")
    return(out)
}



.pema_prep_dat <- function(data, vi, study){
  vi_column <- NULL
  study_column <- NULL
  if(inherits(vi, "character")){
    vi_column <- vi
    vi <- data[[vi]]
    data[[vi_column]] <- NULL
  }
  if(!is.null(study)){
    if(inherits(study, "character")){
      if(!study %in% names(data)) stop("Argument 'study' is not a column of 'data'.")
      study_column <- study
      cl[["study"]] <- data[[study]]
      data[[study_column]] <- NULL
    } else {
      if(!length(study) == nrow(data)){
        stop("Argument 'study' must be a character string referencing a column in 'data', or a vector of study IDs with length equal to the number of rows in 'data'.")
      }
    }
  }
  return(list(vi = vi, study = study, data = data))
}


#' Convert an object to cv.glmnet
#'
#' Create a `cv.glmnet` object from an object for which a method exists,
#' so that all methods for `cv.glmnet` objects can be used.
#' @param x An object for which a method exists.
#' @param ... Arguments passed to or from other methods.
#' @return An object of class `cv.glmnet`,
#' as documented in [glmnet::cv.glmnet()].
#' @export
#' @examples
#' x <- "a"
#' converted <- as.cv.glmnet(x)
as.cv.glmnet <- function(x, ...){
  UseMethod("as.cv.glmnet", x)
}

#' @method as.cv.glmnet pma
#' @export
as.cv.glmnet.pma <- function(x, ...){
  out <- x[["lasso"]]
  class(out) <- c(class(out), "pma_cv_glmnet")
  attr(out, "tau2_cv") <- x$tau2_cv
  attr(out, "tau2") <- x$rma$tau2
  return(out)
}

#' @method as.cv.glmnet character
#' @export
as.cv.glmnet.character <- function(x, ...){
  class(x) <- c("cv.glmnet", "character")
  return(x)
}

#' @method predict pma
#' @export
predict.pma <- function(object, newdata, s = "lambda.min", ...){
  cl <- match.call()
  cl[[1L]] <- quote(predict)
  cl[["object"]] <- as.cv.glmnet(object)
  names(cl)[names(cl) == "newdata"] <- "newx"
  eval.parent(cl)
}

#' @method predict pma
#' @export
predict.pma <- function(object, newdata = NULL, s = "lambda.min", ...){
  cl <- match.call()
  cl[[1L]] <- quote(predict)
  class(object) <- "cv.glmnet"
  cl[["object"]] <- object
  cl[["s"]] <- s
  if(is.null(newdata)){
    cl[["newdata"]] <- object$x
  }
  names(cl)[names(cl) == "newdata"] <- "newx"
  eval.parent(cl)
}


#' @method coef pma
#' @export
coef.pma <- function(object, s = "lambda.min", ...){
  cl <- match.call()
  cl[[1L]] <- quote(coef)
  class(object) <- "cv.glmnet"
  cl[["object"]] <- object
  cl[["s"]] <- s
  eval.parent(cl)
}
