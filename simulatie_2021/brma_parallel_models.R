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

#upper part of brma for simulation function. Note the extra 'method' argument to specify either horseshoe or Lasso
brma_sim <- function(cl, simulated_data, file_stem, method, ...){
  brma_models <- parLapply(cl, simulated_data, function(data) {
    brma(yi ~ .,
         data = data$training,
         method = method,
         ...)
  })
}
