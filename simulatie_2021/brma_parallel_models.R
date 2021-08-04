#load('./simulatie_2021/brma_models.RData') #this loads in brma_models as created in line 44
save(brma_models, file = 'brma_models.RData')

library(parallel)
library(pema) #for brma and simulate_smd
library(brms)  #for package Rcpp that interfaces C++ code into R


#create clusters
cores <- as.integer(floor(detectCores() / 2))
cl <- makeCluster(cores)

#evaluate necessary packages
clusterEvalQ(cl, library(pema))
clusterEvalQ(cl, library(brms))

#create list of simulated datasets to test parallel functions
simulated_data <- list(pema::simulate_smd(), pema::simulate_smd())

#------brma_for_sim function
brma_for_sim <- function(...){
  args <- as.list(match.call()[-1])
  res <- try(do.call(brma, args), silent = TRUE)
}

#load in other necessary functions
source('./simulatie_2021/functions for simulation/model_accuracy.R')

#export functions to clusters
clusterExport(cl, "brma_for_sim")
clusterExport(cl, "model_accuracy")
clusterExport(cl, "pred_brma")

#create brma models. works fine, output is a list with brma() objects
brma_models <- parLapply(cl, simulated_data, function(data) {
  args <- list(formula = force(as.formula(paste0("yi~", paste0(names(data$training)[-c(1:2)], collapse = "+")))),
               data = data$training,
               method = 'lasso')
  do.call(brma_for_sim, args)
})

#also works, returns Boolean matrix of relevant variables for both 95% CI and HDI
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


#does not work....
model_accuracy(
  fit = fit.lasso,
  newdata = as.matrix(simdat$training[, -c(1:2)]),
  observed = simdat[[1]]$yi
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



