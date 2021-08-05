rma_sim <- function(cl, simulated_data, file_stem, ...){
  rma_models <- parLapply(cl, simulated_data, function(data) {
    args <- list(yi = data$training$yi,
                 vi = data$training$vi,
                 mods = as.matrix(data$training[,-c(1:2)]))

    do.call(rma_for_sim, args)
  })
  browser()
  rma_fits <- t(
    clusterMap(
      cl,
      fun = function(models, data) {
        if(is.na(models)){
          return(rep(NA, 7))
        }
        c(
          model_accuracy(
            fit = models,
            newdata = as.matrix(data$training[, -c(1:2)]),
            observed = data[[1]]$yi
          ),
          model_accuracy(
            fit = models,
            newdata = as.matrix(data$testing[, -1]),
            observed = data$testing$yi,
            ymean = mean(data$training$yi)
          ),
          tau2 = models$tau2
        )
      },
      models = rma_models,
      data = simulated_data,
      SIMPLIFY = TRUE,
      USE.NAMES = FALSE
    )
  )

  rma_selected <-
    t(parSapply(
      cl = cl,
      rma_models,
      FUN = function(models) {
        models[["beta"]][-1,1]
      }
    ))

  file_name <-
    paste0(paste(c(file_stem, "fits", chunk), collapse = "_"), ".RData")
  saveRDS(rma_fits, file = file_name, compress = FALSE)

  file_name <-
    paste0(paste(c(file_stem, "selected", chunk), collapse = "_"), ".RData")
  saveRDS(rma_selected,
          file = file_name,
          compress = FALSE)

  rm(rma_models, rma_fits, rma_selected)
}
