#load in libraries
library(pema)
library(brms)

#simulate data
set.seed(6164900)
df <- simulate_smd()
fit <- brma(yi~., data = df$training, method = "lasso")

#following 3 lines return boolean of 95% CI's that do not contain 0
CI95 <- summary(fit$fit)$summary
CI95 <- CI95[paste0("betas[", 1:(ncol(fit$X)-1), "]"), c("2.5%", "97.5%")]
selected <- apply(CI95, 1, function(x){ !(sum(sign(x)) == 0)})

#highest posterior density
hpdi <- bayestestR::hdi(fit$fit, parameters = "betas")
hpdi <- cbind(hpdi$CI_low, hpdi$CI_high)
selected <- apply(CI95, 1, function(x){ !(sum(sign(x)) == 0)})

#load in clusters
cores <- as.integer(floor(detectCores() / 2))
cl <- makeCluster(cores)
clusterEvalQ(cl, library(pema))
clusterEvalQ(cl, library(rstan))
clusterEvalQ(cl, library(brms))
clusterEvalQ(cl, library(sn))

#brma_selected returns matrix of booleans, where True is assigned to relevant predictors
#assessed by both 95%CI (first 5 columns) and HDI (last 5 columns)
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




