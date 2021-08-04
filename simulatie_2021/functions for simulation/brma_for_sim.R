
brma_for_sim <- function(...){
  args <- as.list(match.call()[-1])
  res <- try(do.call(brma, args), silent = TRUE)
}


brma_sim <- function(cl, simulated_data = simulated_data, file_stem, method, ...){
  brma_models <- parLapply(cl, simulated_data, function(data) {
    args <- list(formula = force(as.formula(paste0("yi~", paste0(names(data$training)[-c(1:2)], collapse = "+")))),
                 data = data$training,
                 method = method)
    do.call(brma_for_sim, args)
  })

  #brma fit
  brma_fits <- t(
    clusterMap(
      cl,
      fun = function(models, data) {
        if(is.na(models)){
          return(rep(NA, 7))
        }
        c(
          model_accuracy(
            fit = models,
            newdata = data$training[, -c(1:2)],
            observed = data[[1]]$yi
          ),
          model_accuracy(
            fit = models,
            newdata = data$testing[, -1],
            observed = data$testing$yi,
            ymean = mean(data$training$yi)
          ),
          tau2 = summary(models$fit)$summary['sd_1[1]', 'mean']
        )
      },
      models = brma_models,
      data = simulated_data,
      SIMPLIFY = TRUE,
      USE.NAMES = FALSE
    )
  )

  #brma selected
  brma_selected <- t(
    parSapply(
      cl = cl,
      brma_models,
      FUN = function(models){
        CI95 <- summary(models$fit)$summary[paste0("betas[", 1:(ncol(models$X)-1), "]"), c("2.5%", "97.5%")]
        selected_CI95 <- apply(CI95, 1, function(x){ !(sum(sign(x)) == 0)})
        hpdi <- bayestestR::hdi(models$fit, parameters = "betas")
        hpdi <- cbind(hpdi$CI_low, hpdi$CI_high)
        selected_hpdi <- apply(CI95, 1, function(x){ !(sum(sign(x)) == 0)})
        selected <- cbind(selected_CI95, selected_hpdi)
      }
    )
  )

  file_name <-
    paste0(paste(c(file_stem, "fits", chunk), collapse = "_"), ".RData")
  saveRDS(brma_fits, file = file_name, compress = FALSE)

  file_name <-
    paste0(paste(c(file_stem, "selected", chunk), collapse = "_"), ".RData")
  saveRDS(brma_selected,
          file = file_name,
          compress = FALSE)

  rm(brma_models, brma_fits, brma_selected)
}
