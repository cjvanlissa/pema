#load in libraries
library(pema)
library(brms)

#simulate data
df <- simulate_smd()
fit <- brma(yi~., data = df$training, method = "lasso")

#following 3 lines return boolean of 95% CI's that do not contain 0
CI95 <- summary(fit$fit)$summary
CI95 <- CI95[paste0("betas[", 1:(ncol(fit$X)-1), "]"), c("2.5%", "97.5%")]
selected <- apply(CI95, 1, function(x){ !(sum(sign(x)) == 0)})

#returns 95 CI of relevant beta's, returns empy table if none of variables are relevant
selected_true <- CI95[names(which(selected == T)),]

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

#brma_selected returns EAP of beta_values (just like rma_selected does)
brma_selected <- t(
  parSapply(
    cl = cl,
    brma_models,
    FUN = function(models){
      summary(models$fit)$summary[paste0("betas[", 1:(ncol(models$X)-1), "]"), c("50%")]
    }
  )
)




