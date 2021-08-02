library(pema)
library(brms)
df <- simulate_smd()
fit <- brma(yi~., data = df$training, method = "lasso")

CI95 <- summary(fit$fit)$summary
CI95 <- CI95[paste0("betas[", 1:(ncol(fit$X)-1), "]"), c("2.5%", "97.5%")]
selected <- apply(CI95, 1, function(x){ !(sum(sign(x)) == 0)})

hpdi <- bayestestR::hdi(fit$fit, parameters = "betas")
hpdi <- cbind(hpdi$CI_low, hpdi$CI_high)
selected <- apply(CI95, 1, function(x){ !(sum(sign(x)) == 0)})
