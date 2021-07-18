#' @title Lasso meta-analysis
#' @description Conducts lasso (L1 norm) penalized meta-regression by means of
#' maximum likelihood estimation.
#' @param x Numeric matrix of moderators.
#' @param y Numeric vector of effect sizes.
#' @param v Numeric vector of sampling variances.
#' @param whichweights Character, indicating which meta-analytic weights
#' to use. Default: 'random'
#' @param lambda_n Numeric, indicating how many values of lambda to try.
#' Default: 100
#' @param ... Additional parameters passed to \link[stats]{optim}.
#' @return An object of class 'lma_ml'
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname lma_ml
#' @export
if(FALSE){
  d <- dat.bangertdrowns2004[, c("wic", "feedback", "info", "pers", "imag", "meta", "subject", "yi", "vi")]
  options(na.action='na.pass')


  # Real work starts here ---------------------------------------------------

  # This is a matrix of predictors
  x <- model.matrix(yi~., d)[, -c(25)][, 1:5]
  x[is.na(x)] <- 0
  # Vector of outcomes
  y <- dat.bangertdrowns2004$yi

  # Vector of sampling variances
  v <- dat.bangertdrowns2004$vi
  lma_ml(x, y, v)
  fitGLM <- glmnet(x, y)
  fitGLM$lambda

  library(metaforest)
  library(glmnet)
  d <- SimulateSMD(k_train = 50)
  x <- model.matrix(yi ~ X1+X2+X3+X4+X5, d$training)
  y <- d$training$yi
  v <- d$training$vi
  lma_ml(x, y, v)
  result_glm <- cv.glmnet(x, y)
  result_glm$lambda.min
  coef(result_glm, s = result_glm$lambda.min)
  result_lma <- lma_ml(cbind(1, x), y, v)
  plot(log(1:100), rowMeans(result_lma$cv_results))
  plot(result_glm)
  lambda_n = 100
  standardize = TRUE
  intercept = TRUE
  n_folds = 10
}

lma_ml <- function(x, y, v, whichweights = "random", lambda_n = 100,
                   standardize = TRUE, intercept = TRUE, n_folds = 10, ...){
# Add intercept somewhere
# Check that input is legal -----------------------------------------------
  if(!is.matrix(x)) stop("Argument 'x' must be a matrix.")
  if(!is.numeric(y)) stop("Argument 'y' must be a numeric vector.")
  if(!is.numeric(v)) stop("Argument 'v' must be a numeric vector.")

# Capture input -----------------------------------------------------------
  args_input <- match.call()[-1]

# Standardize predictors --------------------------------------------------
  sd_x <- apply(x, 2, sd_n)
  no_variance <- sd_x == 0
  if(sum(no_variance) > 1){
    stop("Several moderators had no variance.")
  } else {
    if(intercept == TRUE){
      if(sum(no_variance) == 1){
        sd_x[which(no_variance)] <- 1
      } else {
        warning("No intercept found in 'x'. An intercept was added.")
        x <- cbind(1, x)
        sd_x <- c(1, sd_x)
      }
    }
  }

  if(standardize){
    x[, !no_variance] <- scale(x[, !no_variance], scale = sd_x[!no_variance])
  } else {
    sd_x <- rep(1, num_param)
  }

  N <- dim(x)[1]
  num_param <- dim(x)[2]

# Compute range of lambda values --------------------------------------------

  lambda_max <- max(abs(t(x)%*%y))/N

  #lambda_values <- exp(log(10)*seq(lambda_max, 0, length.out = lambda_n))
  epsilon <- .0001
  lambda_values <- round(exp(seq(log(lambda_max), log(lambda_max*epsilon),
                              length.out = lambda_n)), digits = 10)

# Set cv folds ------------------------------------------------------------
  n_folds <- 10
  folds <- sample.int(n_folds, N, replace = TRUE)

# Set optim arguments -----------------------------------------------------
  args_optim <- list(
    par = c(rep(0, num_param), 1e-10),
    fn = lma_loss,
    hessian = TRUE,
    method = "L-BFGS-B",
    lower = c(rep(-Inf, num_param), 1e-10),
    upper = rep(Inf, num_param+1),
    x = x,
    y = y,
    lambda = 0
  )

# Cross-validate analysis -------------------------------------------------
  cv_results <- sapply(1:n_folds, function(k) {
    fold <- which(folds == k)
    #fold <- which(folds == 1)
    args_optim$x <- x[-fold,]
    args_optim$y <- y[-fold]
    sapply(lambda_values, function(this_lambda) {
      #this_lambda <- lambda_values[4]
      args_optim$lambda <- this_lambda

      this_fit <- do.call(optim, args_optim)
      resid <- y[fold] - (x[fold, ] %*% head(this_fit$par,-1))
      mean(resid ^ 2)
    })
  })


  args_optim$lambda <- lambda_values[which.min(rowMeans(cv_results))]
  parameter.fits <- do.call(optim, args_optim)
  browser()
  results <- lma_lm_results_table(parameter.fits, sd_x)
  #rownames(results) <- c(colnames(x), "tau2")
  out <- list(results = results,
              tau2 = parameter.fits$par[length(parameter.fits$par)],
              lambda.min = lambda_values[which.min(rowMeans(cv_results))],
              cv_results = cv_results,
              lambda_values = lambda_values
              )
  class(out) <- "lma_ml"
  out
}

# Loss function for lasso penalized random effects model
lma_loss <- function(par, x, y, lambda) {
  beta <- par[-length(par)]
  tau2 <- par[length(par)]
  #if(tau2 == 0) return(1e12)
  resid <- y - (x %*% beta)

  log.likelihoods <- dnorm(resid, mean = 0, sd = sqrt(tau2), log = TRUE)
  -sum(log.likelihoods) + (lambda * sum(abs(beta)))
}

# Function to get SD dividing by N instead of N-1
sd_n <- function(x) sqrt(sum((x-mean(x))^2)/length(x))

# Function to get nicely formatted results table
lma_lm_results_table <- function(x, sds = NULL){
  hessian.inv <- solve(x$hessian)
  parameter.se <- sqrt(diag(hessian.inv))
  Z = x$par/parameter.se

  results <- cbind(est = x$par,
                   se = parameter.se,
                   Z = Z,
                   p = 2*pnorm(abs(Z), lower.tail = FALSE),
                   ci_lo = x$par - 1.96 * parameter.se,
                   ci_hi = x$par + 1.96 * parameter.se)
}
