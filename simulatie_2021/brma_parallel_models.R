library(parallel)
library(metaforest) #for random forest
library(pema) #for brma and simulate_smd
library(rstan) #for compiling C++
library(brms)  #for package Rcpp that interfaces C++ code into R
library(sn)

#create clusters
cores <- as.integer(floor(detectCores() / 2))
cl <- makeCluster(cores)

#evaluate necessary packages
clusterEvalQ(cl, library(metaforest))
clusterEvalQ(cl, library(pema))
clusterEvalQ(cl, library(rstan))
clusterEvalQ(cl, library(brms))
clusterEvalQ(cl, library(sn))


#create list of simulated datasets to test parallel functions
simulated_data <- list(pema::simulate_smd(), pema::simulate_smd())

#first test on metaforest with 100 trees. works.
mf_r_models <- parLapply(cl, simulated_data, function(data) {
  MetaForest(yi ~.,
              data = data$training,
              num.trees = 100)
  })

#test on lasso model. Works! (but only still if variables "vi" and "study" are explicitly defined)
lasso_models <- parLapply(cl, simulated_data, function(data) {
  brma(yi ~ .,
       data = data$training,
       vi = data$training$vi,
       study = 1:nrow(data$training),
       method = 'lasso')})

#test on horseshoe model. Works! (but only still if variables "vi" and "study" are explicitly defined)
hs_models <- parLapply(cl, simulated_data, function(data) {
  brma(yi ~ .,
       data = data$training,
       vi = data$training$vi,
       study = 1:nrow(data$training),
       method = 'hs')})


#stop clusters
stopCluster(cl)


#------work from here
brma_for_sim <- function(...){
  args <- as.list(match.call()[-1])
  res <- try(do.call(brma, args), silent = TRUE)
}

clusterExport(cl, "brma_for_sim")

#upper part of brma for simulation function. Note the extra 'method' argument to specify either horseshoe or Lasso
#this part works
brma_sim <- function(cl, simulated_data, file_stem, method, ...){
  brma_models <- parLapply(cl, simulated_data, function(data) {
    args <- list(formula = force(as.formula(paste0("yi~", paste0(names(data$training)[-c(1:2)], collapse = "+")))),
                            vi = data$training$vi,
                            study = 1:nrow(data$training),
                            data = data$training,
                            method = 'lasso')
    do.call(brma_for_sim, args)
  })
}


rma_sim <- function(cl, simulated_data, file_stem, ...){
  rma_models <- parLapply(cl, simulated_data, function(data) {
    args <- list(yi = data$training$yi,
                 vi = data$training$vi,
                 mods = as.matrix(data$training[,-c(1:2)]))

    do.call(rma_for_sim, args)
  })

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
