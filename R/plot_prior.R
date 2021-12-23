plot_prior <- function(method = c("hs", "lasso"),
                       prior = switch(method,
                                      "lasso" = c(df = 1, scale = 1),
                                      "hs" = c(df = 1, df_global = 1, df_slab = 4, scale_global = 1, scale_slab = 1, par_ratio = NULL)),
                       iter = 1000){
  smpls <- suppressWarnings(sampling(object = pema:::stanmodels[[c("lasso_prior", "hs_prior")[(method[1] == "hs")+1]]], data = as.list(prior), chains = 1, iter = 10000, warmup = 0, show_messages = FALSE, verbose = FALSE, refresh = 0))
  plot(density(smpls@sim$samples[[1]]$b), main = c("Lasso prior", "Horseshoe prior")[(method == "hs")+1],
       xlab = paste0("Samples: ", iter, ", ", paste0(names(prior), " = ", prior, collapse = ", ")),
       xlim = c(-5, 5))
}

## Plotting the shrinkage priors in pema
library(ggplot2)
library(tidyr)
library(rmutil)
library(LaplacesDemon)

set.seed(23122021)

ndraws <- 1e+05

## lasso
lasso <- rmutil::rlaplace(ndraws, m=0, s=0.5)

## regularized horseshoe
reg.hs <- rep(NA, ndraws)
for(i in 1:ndraws){
  c2 <- LaplacesDemon::rinvgamma(1, shape=0.5, scale=1)
  lambda <- LaplacesDemon::rhalfcauchy(1, scale=1)
  tau <- LaplacesDemon::rhalfcauchy(1, scale=1)
  lambda2_tilde <- c2 * lambda^2/(c2 + tau^2*lambda^2)
  reg.hs[i] <- rnorm(1, 0, sqrt(tau^2*lambda2_tilde))
}

## plot
df.comb <- data.frame(lasso, reg.hs)
df.long <- gather(df.comb, Prior, value) # long format
df.long$Prior <- factor(df.long$Prior)
levels(df.long$Prior) <- list("Lasso"="lasso", "Horseshoe"="reg.hs")

p1 <- ggplot(df.long, aes(x = value, linetype = Prior)) +
  stat_density(geom = "line", position = "identity") +
  scale_linetype_manual(values = c(1, 2))
p1

# xlim removes values outside the range, coord_cartesian does not, but behaves strange sometimes
# solution: restrict range with xlim first & then use coord-cartesian
p2 <- ggplot(df.long, aes(x = value, linetype = Prior)) +
  stat_density(geom = "line", position = "identity") +
  scale_linetype_manual(values = c(1, 2)) +
  xlim(-10,10) + coord_cartesian(xlim=c(-5,5))
p2
