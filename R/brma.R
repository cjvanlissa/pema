#' Conduct Bayesian Regularized Meta-Analysis
#'
#' This function conducts Bayesian regularized meta-regression (Van Lissa & Van
#' Erp, 2021). It uses the \code{stan} function
#' [rstan::sampling] to fit the model. A lasso or horseshoe prior is used to
#' shrink the regression coefficients of irrelevant moderators towards zero.
#' See Details.
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
#' @param standardize Logical, indicating whether or not to standardize the
#' predictors (defaults to \code{TRUE}, which is recommended so that shrinking
#' affects all parameters similarly.
#' @param prior Numeric vector, specifying the prior to use. Note that the
#' different \code{method}s require this vector to contain specific named
#' elements.
# @param iter A positive integer specifying the number of iterations for each
# chain (including warmup). Defaults to 2000.
#  the model statement. Defaults to .5.
# @param chains A positive integer specifying the number of Markov chains.
# Defaults to 4.
#' @param mute_stan Logical, indicating whether mute all 'Stan' output or not.
#' @param ... Additional arguments passed on to [rstan::sampling()].
#' Use this, e.g., to override default arguments of that function.
#' @details The Bayesian regularized meta-analysis algorithm (Van Lissa & Van
#' Erp, 2021) penalizes meta-regression coefficients either via the
#' lasso prior (Park & Casella, 2008) or the regularized horseshoe prior
#' (Piironen & Vehtari, 2017).
#' \describe{
#'   \item{lasso}{ The Bayesian equivalent of the lasso penalty is obtained when
#'   placing independent Laplace (i.e., double exponential) priors on the
#'   regression coefficients centered around zero. The scale of the Laplace
#'   priors is determined by a global scale parameter \code{scale}, which
#'   defaults to 1 and an inverse-tuning parameter \eqn{\frac{1}{\lambda}}
#'   which is given a chi-square prior governed by a degrees of freedom
#'   parameter \code{df} (defaults to 1). If \code{standardize = TRUE},
#'   shrinkage will
#'   affect all coefficients equally and it is not necessary to adapt the
#'   \code{scale} parameter. Increasing the \code{df} parameter will allow
#'   larger values for the inverse-tuning parameter, leading to less shrinkage.}
#'   \item{hs}{ One issue with the lasso prior is that it has relatively light
#'   tails. As a result, not only does the lasso have the desirable behavior of
#'   pulling small coefficients to zero, it also results in too much shrinkage
#'   of large coefficients. An alternative prior that improves upon this
#'   shrinkage pattern is the horseshoe prior (Carvalho, Polson & Scott, 2010).
#'   The horseshoe prior has an infinitely large spike at zero, thereby pulling
#'   small coefficients toward zero but in addition has fat tails, which allow
#'   substantial coefficients to escape the shrinkage. The regularized horseshoe
#'   is an extension of the horseshoe prior that allows the inclusion of prior
#'   information regarding the number of relevant predictors and can
#'   be more numerically stable in certain cases (Piironen & Vehtari, 2017).
#'   The regularized horseshoe has a global shrinkage parameter that influences
#'   all coefficients similarly and local shrinkage parameters that enable
#'   flexible shrinkage patterns for each coefficient separately. The local
#'   shrinkage parameters are given a Student's t prior with a default \code{df}
#'   parameter of 1. Larger values for \code{df} result in lighter tails and
#'   a prior that is no longer strictly a horseshoe prior. However, increasing
#'   \code{df} slightly might be necessary to avoid divergent transitions in
#'   Stan (see also \url{https://mc-stan.org/misc/warnings.html}). Similarly,
#'   the degrees of freedom for the Student's t prior on the global shrinkage
#'   parameter \code{df_global} can be increased from the default of 1 to, for
#'   example, 3 if divergent transitions occur although the resulting
#'   prior is then strictly no longer a horseshoe. The scale for the Student's t
#'   prior on the global shrinkage parameter \code{scale_global} defaults to 1
#'   and can be decreased to achieve more shrinkage. Moreover, if prior
#'   information regarding the number of relevant moderators is available, it is
#'   recommended to include this information via the \code{par_ratio} argument
#'   by setting it to the ratio of the expected number of non-zero coefficients
#'   to the expected number of zero coefficients. When \code{par_ratio} is
#'   specified, \code{scale_global} is ignored and instead based on the
#'   available prior information. Contrary to the horseshoe prior, the
#'   regularized horseshoe applies additional regularization on large
#'   coefficients which is governed by a Student's t prior with a
#'   \code{scale_slab} defaulting to 2 and \code{df_slab} defaulting to 4.
#'   This additional regularization ensures at least some shrinkage of large
#'   coefficients to avoid any sampling problems.}
#' }
#' @references
#' Van Lissa, C. J., & van Erp, S. (2021, December 9). Select relevant
#' moderators using Bayesian regularized meta-regression.
#' \doi{10.31234/osf.io/6phs5}

#' Park, T., & Casella, G. (2008). The Bayesian Lasso. Journal of the American
#' Statistical Association, 103(482), 681–686. \doi{10.1198/016214508000000337}
#'
#' Carvalho, C. M., Polson, N. G., & Scott, J. G. (2010). The horseshoe
#' estimator for sparse signals. Biometrika, 97(2), 465–480.
#' \doi{10.1093/biomet/asq017}
#'
#' Piironen, J., & Vehtari, A. (2017). Sparsity information and regularization
#' in the horseshoe and other shrinkage priors. Electronic Journal of
#' Statistics, 11(2). \url{https://projecteuclid.org/journals/electronic-journal-of-statistics/volume-11/issue-2/Sparsity-information-and-regularization-in-the-horseshoe-and-other-shrinkage/10.1214/17-EJS1337SI.pdf}
#' @return A `list` object of class `brma`, with the following structure:
#' ```
#' list(
#'   fit          # An object of class stanfit, for compatibility with rstan
#'   coefficients # A numeric matrix with parameter estimates; these are
#'                # interpreted as regression coefficients, except tau2 and tau,
#'                # which are interpreted as the residual variance and standard
#'                # deviation, respectively.
#'   formula      # The formula used to estimate the model
#'   terms        # The predictor terms in the formula
#'   X            # Numeric matrix of moderator variables
#'   Y            # Numeric vector with effect sizes
#'   vi           # Numeric vector with effect size variances
#'   tau2         # Numeric, estimated tau2
#'   R2           # Numeric, estimated heterogeneity explained by the moderators
#'   k            # Numeric, number of effect sizes
#'   vi_column    # Optional, name of the column in the original data
#'                # corresponding to vi
#'   study        # Numeric vector with study id numbers
#'   study_column # Optional, name of the column in the original data
#'                # corresponding to study
#' )
#' ```
#' @export
#' @examples
#' data("curry")
#' df <- curry[c(1:5, 50:55), c("d", "vi", "sex", "age", "donorcode")]
#' suppressWarnings({res <- brma(d~., data = df, iter = 10)})
#' @importMethodsFrom rstan summary
#' @importFrom stats model.matrix na.omit quantile sd
#' @importFrom RcppParallel RcppParallelLibs CxxFlags
#' @importFrom rstantools bayes_R2
# The line above is just to avoid CRAN warnings that RcppParallel is not
# imported from, despite RcppParallel being a necessary dependency of rstan.
brma <-
  function(formula,
           data,
           vi = "vi",
           study = NULL,
           method = "hs",
           standardize = TRUE,
           prior = switch(method,
                          "lasso" = c(df = 1, scale = 1),
                          "hs" = c(df = 1, df_global = 1, df_slab = 4, scale_global = 1, scale_slab = 1, par_ratio = NULL)),
           mute_stan = TRUE,
           ...) {
    # Bookkeeping for columns that should not be in X or Y
    vi_column <- NULL
    study_column <- NULL
    if(inherits(vi, "character")){
      vi_column <- vi
      vi <- data[[vi]]
      data[[vi_column]] <- NULL
    }
    if(is.null(study)){
      study <- 1:nrow(data)
    } else {
      if(inherits(study, "character")){
        study_column <- study
        study <- data[[study]]
        data[[study_column]] <- NULL
      }
    }
    # Make model matrix
    mf <- match.call(expand.dots = FALSE)
    mf <- mf[c(1L, match(c("formula", "subset", "na.action"), names(mf), 0L))]
    mf[["data"]] <- data
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- quote(stats::model.frame)
    mf <- eval(mf, parent.frame())
    Y <- mf[[1]]
    #X <- mf[, -1, drop = FALSE]
    mt <- attr(mf, "terms")
    X <- model.matrix(mt, mf)
    se <- sqrt(vi)
    N <- length(Y)
    # if(isTRUE(standardize)){
    #   X <- scale(X) # Should there be any fancy standardization for categorical variables?
    #                 # Should coefficients be transformed back to original scale?
    #   scale_m <- attr(X, "scaled:center")
    #   scale_s <- attr(X, "scaled:scale")
    # }
    #X <- cbind(1, X)
    standat <- c(
      list(
        N = N,
        Y = Y,
        se = se,
        K = ncol(X),
        X = X),
      as.list(prior),
      list(
        N_1 = length(unique(study)),
        M_1 = 1,
        J_1 = study,
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
    # Mute stan
    dots <- list(...)
    if(!any(c("show_messages", "verbose", "refresh") %in% names(dots))){
      if(mute_stan){
        cl[["show_messages"]] <- FALSE
        cl[["verbose"]] <- FALSE
        cl[["refresh"]] <- 0
      }
    }
    fit <- eval(cl)
    sums <- summary(fit)$summary
    keepthese <- c(which(rownames(sums) == "Intercept"),
                   which(startsWith(rownames(sums), "b[")),
                   which(rownames(sums) == "sd_1[1]"))
    sums <- sums[keepthese, , drop = FALSE]
    sim <- fit@sim
    sdpar <- match("sd_1[1]", sim$fnames_oi)
    tau2 <- .extract_samples(sim, sdpar)^2
    addrow <- sums["sd_1[1]",]
    addrow[c("mean", "sd", "2.5%", "25%", "50%", "75%", "97.5%")] <-
      c(mean(tau2), sd(tau2), quantile(tau2, c(.025, .25, .5, .75,.975)))
    addrow["se_mean"] <- addrow["sd"]/sqrt(addrow[["n_eff"]])
    #tau <- unlist(lapply(fit@sim$samples, `[[`, "sd_1[1]"))
    sums <- rbind(sums, tau2 = addrow)
    rownames(sums)[2:(ncol(X))] <- colnames(X)[-1]
    rownames(sums)[rownames(sums) == "sd_1[1]"] <- "tau"
    tau2 <- unname(addrow["mean"])
    Wi <- 1 / vi
    tau2_before <-
      max(0, (sum(Wi * (Y - (
        sum(Wi * Y) / sum(Wi)
      )) ^ 2) - (N - 1)) / (sum(Wi) - (sum(Wi ^ 2) / sum(Wi))))

    R2 <- max(0, 100 * (tau2_before-tau2)/tau2_before)
    #I2 <- 100 * tau2/(vt + tau2)
    #H2 <- tau2/vt + 1

    out <- list(fit = fit,
                coefficients = sums,
                formula = formula,
                terms = mt,
                X = X,
                Y = Y,
                vi = vi,
                tau2 = tau2,
                #I2 = I2,
                #H2 = H2,
                R2 = R2,
                k = N)
    if(!is.null(vi_column)) out$vi_column <- vi_column
    if(!is.null(study_column)) out$study_column <- study_column
    if(!is.null(study)) out$study <- study
    class(out) <- c("brma", class(out))
    return(out)
  }
