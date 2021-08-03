#load('./simulatie_2021/brma_models.RData') #this loads in brma_models as created in line 44

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


#------brma_for_sim function
brma_for_sim <- function(...){
  args <- as.list(match.call()[-1])
  res <- try(do.call(brma, args), silent = TRUE)
}

#load in other necessary functions
source('./simulatie_2021/functions for simulation/model_accuracy.R')
source('./simulatie_2021/functions for simulation/rma_sim.R')
source('./simulatie_2021/functions for simulation/rma_for_sim.R')

#export functions to clusters
clusterExport(cl, "brma_for_sim")
clusterExport(cl, "model_accuracy")
clusterExport(cl, "pred_brma")
clusterExport(cl, "rma_for_sim")
clusterExport(cl, "rma_sim")

#create brma models. works fine, output is a list with brma() objects
brma_models <- parLapply(cl, simulated_data, function(data) {
  args <- list(formula = force(as.formula(paste0("yi~", paste0(names(data$training)[-c(1:2)], collapse = "+")))),
               data = data$training,
               method = 'lasso')
  do.call(brma_for_sim, args)
})

#also works, returns EAP of beta values (just like output for rma_selected)
brma_selected <- t(
  parSapply(
    cl = cl,
    brma_models,
    FUN = function(models){
      summary(models$fit)$summary[paste0("betas[", 1:(ncol(models$X)-1), "]"), c("50%")]
    }
  )
)

#creates error: 'replacement has length zero'
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
    models = brma_models,
    data = simulated_data,
    SIMPLIFY = TRUE,
    USE.NAMES = FALSE
  )
)





### -------- code below is ongoing proces, not necessary to look at.

#works up until do.call on line 46. Note the extra 'method' argument to specify either horseshoe or Lasso
brma_sim <- function(cl, simulated_data, file_stem, method, ...){
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
      models = brma_models,
      data = simulated_data,
      SIMPLIFY = TRUE,
      USE.NAMES = FALSE
    )
  )
}

#stop clusters
stopCluster(cl)

out2 <- pred_brma(brma_models[1])

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
