run_rma <- function(data, max_mods = 6) {
  out <- rep(NA, max_mods)
  tstart <- Sys.time()
  models <- tryCatch(
    metafor::rma(
      yi = data$training$yi,
      vi = data$training$vi,
      mods = as.matrix(data$training[, -c(1:2)])
    ),
    error = function(e) NULL
  )
  runtime <- as.numeric(Sys.time() - tstart)
  if (!is.null(models)) {
    fits <- c(
      model_accuracy(
        fit = models,
        newdata = as.matrix(data$training[, -c(1:2)]),
        observed = data$training$yi
      ),
      model_accuracy(
        fit = models,
        newdata = as.matrix(data$testing[, -1]),
        observed = data$testing$yi,
        ymean = mean(data$training$yi)
      ),
      tau2 = models$tau2,
      beta1 = models$b[2, 1]
    )
    out[seq_along(models$pval[-1])] <- models$pval[-1] < 0.05

  } else {
    fits <- rep(NA, 4)
  }
  c(fits, out, runtime)
}

run_2step <- function(data, max_mods = 6) {
  out <- rep(NA, max_mods)
  tstart <- Sys.time()
  models <- tryCatch({
    fit_las <- glmnet::cv.glmnet(
      x = as.matrix(data$training[, -c(1:2)]),
      y = data$training$yi,
      alpha = 1,
      nfolds = 5
    )
  }, error = function(e){NULL})
  runtime <- as.numeric(Sys.time() - tstart)
  if(!is.null(models)){
    mods_min <- as.matrix(
      data$training[
        , -c(
          1, 2,
          which(coef(fit_las, s = "lambda.min")[-1, 1] == 0) + 2
        ),
        drop = FALSE
      ]
    )
    res_min <- suppressWarnings(metafor::rma(
      yi = data$training$yi,
      vi = data$training$vi,
      mods = mods_min
    ))
    mods_1se <- as.matrix(
      data$training[
        , -c(
          1, 2,
          which(coef(fit_las, s = "lambda.1se")[-1, 1] == 0) + 2
        ),
        drop = FALSE
      ]
    )
    res_1se <- suppressWarnings(metafor::rma(
      yi = data$training$yi,
      vi = data$training$vi,
      mods = mods_1se
    ))
    fits_min <- c(
      model_accuracy(
        fit = res_min,
        newdata = res_min$X[, -1, drop = FALSE],
        observed = data$training$yi
      ),
      model_accuracy(
        fit = res_min,
        newdata = as.matrix(
          data$testing[, rownames(res_min$b)[-1], drop = FALSE]
        ),
        observed = data$testing$yi,
        ymean = mean(data$training$yi)
      ),
      tau2 = res_min$tau2,
      beta1 = tryCatch(
        res_min$b["X1", 1],
        error = function(e) 0
      )
    )
    fits_1se <- c(
      model_accuracy(
        fit = res_1se,
        newdata = res_1se$X[, -1, drop = FALSE],
        observed = data$training$yi
      ),
      model_accuracy(
        fit = res_1se,
        newdata = as.matrix(
          data$testing[, rownames(res_1se$b)[-1], drop = FALSE]
        ),
        observed = data$testing$yi,
        ymean = mean(data$training$yi)
      ),
      tau2 = res_1se$tau2,
      beta1 = tryCatch(
        res_1se$b["X1", 1],
        error = function(e) 0
      )
    )

    out_1se <- out_min <- rep(NA, max_mods)
    out_min[1:(ncol(data$training)-2L)] <- names(data$training)[-c(1:2)] %in% rownames(res_min$b)
    out_1se[1:(ncol(data$training)-2L)] <- names(data$training)[-c(1:2)] %in% rownames(res_1se$b)
    out <- c(fits_min, out_min, fits_1se, out_1se, runtime)
  } else {
    out <- rep(NA, sum(c(4, max_mods, 4, max_mods, 1)))
  }
  return(out)
}

run_brma_las <- function(data, max_mods = 6) {

  out <- rep(NA, max_mods)
  tstart <- Sys.time()
  models <- tryCatch(
    brma(
      formula = yi ~ .,
      data = data$training,
      method = "lasso",
      chains = 1,
      iter = 6000,
      refresh = -1,
      mute_stan = TRUE
    ),
    error = function(e) NULL
  )
  runtime <- as.numeric(Sys.time() - tstart)
  if (!is.null(models)) {
    summ <- summary(models$fit)$summary

    fits <- c(
      model_accuracy(
        fit = models,
        newdata = data$training[, -1],
        observed = data$training$yi
      ),
      model_accuracy(
        fit = models,
        newdata = data$testing,
        observed = data$testing$yi,
        ymean = mean(data$training$yi)
      ),
      tau2 = median(models$fit@sim$samples[[1]]$`sd_1[1]`^2),
      beta1 = median(models$fit@sim$samples[[1]]$`betas[1]`),
      divergent = sum(
        rstan:::get_sampler_params(models$fit)[[1]][
          3001:6000, "divergent__", drop = TRUE
        ]
      ),
      maxRhat = max(summ[, "Rhat"], na.rm = TRUE),
      minNeff = min(summ[, "n_eff"], na.rm = TRUE)
    )

    hpdi <- apply(sapply(1:(ncol(data$training)-2L), function(i) pema:::hdi(pema:::.extract_samples(models$fit@sim, paste0("betas[", i, "]")))), 2, function(x) !(sum(sign(x)) == 0))

    out[1:length(hpdi)] <- hpdi


  } else {
    return(rep(NA, 14))
  }
  return(c(fits, out, runtime))
}


run_brma_hs <- function(data, max_mods = 6) {

  out <- rep(NA, max_mods)
  tstart <- Sys.time()
  models <- tryCatch(
    brma(
      formula = yi ~ .,
      data = data$training,
      method = "hs",
      chains = 1,
      iter = 6000,
      refresh = -1,
      verbose = FALSE,
      mute_stan = TRUE
    ),
    error = function(e) NULL
  )
  runtime <- as.numeric(Sys.time() - tstart)
  if (!is.null(models)) {
    summ <- summary(models$fit)$summary

    fits <- c(
      model_accuracy(
        fit = models,
        newdata = data$training[, -1],
        observed = data$training$yi
      ),
      model_accuracy(
        fit = models,
        newdata = data$testing,
        observed = data$testing$yi,
        ymean = mean(data$training$yi)
      ),
      tau2 = median(models$fit@sim$samples[[1]]$`sd_1[1]`^2),
      beta1 = median(models$fit@sim$samples[[1]]$`betas[1]`),
      divergent = sum(
        rstan:::get_sampler_params(models$fit)[[1]][
          3001:6000, "divergent__", drop = TRUE
        ]
      ),
      maxRhat = max(summ[, "Rhat"], na.rm = TRUE),
      minNeff = min(summ[, "n_eff"], na.rm = TRUE)
    )

    hpdi <- apply(sapply(1:(ncol(data$training)-2L), function(i) pema:::hdi(pema:::.extract_samples(models$fit@sim, paste0("betas[", i, "]")))), 2, function(x) !(sum(sign(x)) == 0))

    out[1:length(hpdi)] <- hpdi


  } else {
    return(rep(NA, 18))
  }
  return(c(fits, out, runtime))
}



run_lrma <- function(data, max_mods = 6) {
  out <- rep(NA, max_mods)
  tstart <- Sys.time()
  models <- tryCatch(
    pema:::cv.lrma.formula(
      formula = yi ~ .,
      data = data$training
    ),
    error = function(e) NULL
  )
  runtime <- as.numeric(Sys.time() - tstart)

  if (!is.null(models)) {

    fits_min <- c(
      model_accuracy(
        fit = models,
        newdata = data$training[, -1],
        observed = data$training$yi,
        lambda = "min"
      ),
      model_accuracy(
        fit = models,
        newdata = data$testing,
        observed = data$testing$yi,
        ymean = mean(data$training$yi),
        lambda = "min"
      ),
      tau2 = models$fit$tau2,
      beta1 = models$fit$coefficients[2]
    )
    fits_1se <- c(
      model_accuracy(
        fit = models,
        newdata = data$training[, -1],
        observed = data$training$yi,
        lambda = "1se"
      ),
      model_accuracy(
        fit = models,
        newdata = data$testing,
        observed = data$testing$yi,
        ymean = mean(data$training$yi),
        lambda = "1se"
      ),
      tau2 = models$fit.1se$tau2,
      beta1 = models$fit.1se$coefficients[2]
    )

    out_1se <- out_min <- rep(NA, max_mods)
    out_min[1:(ncol(data$training)-2L)] <- !abs(models$fit$coefficients[-1]) < 1e-4
    out_1se[1:(ncol(data$training)-2L)] <- !abs(models$fit.1se$coefficients[-1]) < 1e-4
    out <- c(fits_min, out_min, fits_1se, out_1se, runtime)
  } else {
    out <- rep(NA, sum(c(4, max_mods, 4, max_mods, 1)))
  }
  return(out)
}


