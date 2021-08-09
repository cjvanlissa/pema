brma_sim <- function(cl, simulated_data = simulated_data, file_stem, method, ...){
  brma_models <- parLapply(cl, simulated_data, function(data) {
    tryCatch({brma(formula = yi~.,
                   data = data$training,
                   method = method,
                   chains = 1,
                   iter = 6000)},
             error = function(e){NULL})
  })
  brma_fits <- t(
    clusterMap(
      cl,
      fun = function(models, data) {
        if(is.null(models)){
          return(rep(NA, 11))
        } else {
          summ <- summary(models$fit)$summary
          c(
            model_accuracy(
              fit = models,
              newdata = data$training[, -1],
              observed = data[[1]]$yi
            ),
            model_accuracy(
              fit = models,
              newdata = data$testing,
              observed = data$testing$yi,
              ymean = mean(data$training$yi)
            ),
            tau2 = mean(models$fit@sim$samples[[1]]$`sd_1[1]`),
            beta1 = mean(models$fit@sim$samples[[1]]$betas[1]),
            divergent = {
              np = nuts_params(models$fit); sum(np$Parameter == "divergent__" & np$Value == 1)
            },
            maxRhat = max(summ[, "Rhat"], na.rm = TRUE),
            minNeff = min(summ[, "n_eff"], na.rm = TRUE)
          )
        }
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
      CI95 <- summary(models$fit)$summary[paste0("betas[", 1:(ncol(models$X)-1), "]"), c("2.5%", "97.5%"), drop = FALSE]
      selected_CI95 <- apply(CI95, 1, function(x){ !(sum(sign(x)) == 0)})
      hpdi <- bayestestR::hdi(models$fit, parameters = "betas")
      hpdi <- cbind(hpdi$CI_low, hpdi$CI_high)
      selected_hpdi <- apply(CI95, 1, function(x){ !(sum(sign(x)) == 0)})
      c(selected_CI95, selected_hpdi)
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
