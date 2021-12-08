df <- simulate_smd(k_train = 500)$training
fit <- brma(yi~., data = df, iter = 100, seed = 30)

CI95 <- summary(fit$fit)$summary
CI95 <- CI95[paste0("betas[", 1:(ncol(fit$X)-1), "]"), c("2.5%", "97.5%")]
selected <- apply(CI95, 1, function(x){ !(sum(sign(x)) == 0)})

test_that("ci selects relevant vars", {
  expect_true(selected[1])
  expect_true(all(!selected[-1]))
})
#hpdi <- bayestestR::hdi(fit$fit, parameters = "betas")
#hpdi <- cbind(hpdi$CI_low, hpdi$CI_high)
#selected <- apply(CI95, 1, function(x){ !(sum(sign(x)) == 0)})
fitrma <- rma(df$yi, df$vi)
res <- rma(object$Y, object$vi, mods = object$X[, -1])
