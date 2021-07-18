##### BAYESIAN REGULARIZED META-ANALYSIS -----
## Author: Sara van Erp
# Some code to run Bayesian regularized meta-analysis with the lasso and regularized horseshoe prior using brms.

library(brms)
library(projpred)
library(bayesplot)
library(sn) # for rsn function

set.seed(280621)

##### Simulate some data using Caspar's function -----
## lma function relies on "rsn" which I assume is to simulate from a skew-normal
## I adapted the function to compute the SE as well for the weighting of studies
simulate_smd <- function(k_train = 20, k_test = 100, mean_n = 40, es = .5,
                         tau2 = 0.04, alpha, moderators = 5, distribution = "normal",
                         model = "es * x[, 1]")
{
  # Make label for training and test datasets
  training <- c(rep(1, k_train), rep(0, k_test))

  # Randomly determine n from a normal distribution with mean
  # mean_n and sd mean_n/3
  n <- rnorm(length(training), mean = mean_n, sd = mean_n/3)
  # Round to even number
  n <- ceiling(n) - ceiling(n)%%2
  # Truncate n so the lowest value can be 8; cases where n=2 will result
  # in errors
  n[n < 8] <- 8

  # Generate moderator matrix x:
  if(distribution == "normal") x <- matrix(rnorm(length(n) * moderators), ncol = moderators)
  if(distribution == "bernoulli") x <- matrix(rbinom((length(n) * moderators), 1, .5), ncol = moderators)
  if(!(distribution %in% c("normal", "bernoulli"))) stop(paste0(distribution, "is not a valid distribution for SimulateSMD"))

  # Sample true effect sizes theta.i from a normal distribution with mean
  # mu, and variance tau2, where mu is the average
  # population effect size. The value of mu depends on the values of the
  # moderators and the true model mu <- eval(model)
  model <- parse(text = model)
  mu <- eval(model)

  # theta.i: true effect size of study i
  theta.i <- mu + rsn(n = length(n), xi = 0, omega = sqrt(tau2), alpha = alpha)

  # Then the observed effect size yi is sampled from a non-central
  # t-distribution under the assumption that the treatment group and
  # control group are both the same size
  p_ntk <- 0.5  #Percentage of cases in the treatment group
  ntk <- p_ntk * n  #n in the treatment group for study i
  nck <- (1 - p_ntk) * n  #n in the control group for study i
  df <- n - 2  #degrees of freedom
  j <- 1 - 3/(4 * df - 1)  #correction for bias
  nk <- (ntk * nck)/(ntk + nck)
  ncp <- theta.i * sqrt(nk)  #Non-centrality parameter

  # Standardized mean difference drawn from a non-central t-distribution
  SMD <- mapply(FUN = rt, n = 1, df = df, ncp = ncp)

  # yi is Hedges' g for study i
  yi <- SMD/((j^-1) * (nk^0.5))

  # Calculate the variance of the effect size
  vi <- j^2 * (((ntk + nck)/(ntk * nck)) + ((yi/j)^2/(2 * (ntk + nck))))

  # Calculate the SE
  SEi <- sqrt(vi/(ntk+nck)) # Why vi divided by ntk+nck and not just sqrt(vi)?

  # Dersimonian and Laird estimate of tau2
  Wi <- 1/vi[1:k_train]
  tau2_est <- max(0, (sum(Wi * (yi[1:k_train] - (sum(Wi * yi[1:k_train])/sum(Wi)))^2) -
                        (k_train - 1))/(sum(Wi) - (sum(Wi^2)/sum(Wi))))

  data <- data.frame(training, vi, SEi, yi, x)

  list(training = subset(data, training == 1, -1), testing = subset(data,
                                                                    training == 0, -c(1, 2)), housekeeping = data.frame(n = n, mu_i = mu, theta_i = theta.i),
       tau2_est = tau2_est)
}

simdat <- simulate_smd(alpha = 0) # using the defaults: 5 moderators with only the 1st being relevant

# add info on the study from which an effect size comes
simdat$training$study <- 1:nrow(simdat$training)

##### Fit the mixed-effects meta-regression model using brms -----
## based on: https://bookdown.org/MathiasHarrer/Doing_Meta_Analysis_in_R/bayesian-meta-analysis-in-r-using-the-brms-package.html

# moderators should be on the same scale to ensure that the shrinkage affects all coefficients equally
# standardize the moderators
cols <- sapply(paste0("X", 1:5), function(x) grep(x, colnames(simdat$training)))
scaled.mods <- scale(simdat$training[, cols])
scaled.training <- cbind.data.frame(scaled.mods, simdat$training[, -cols])

# check the default priors
get_prior(yi|se(SEi) ~ . -vi - study + (1|study), data = scaled.training)
# the intercept gets a student_t prior with df = 3, mean = median(y), and scale 2.5 based on some cutoff for y (see here: https://discourse.mc-stan.org/t/default-student-t-priors-in-brms/17197/8)
# the standard deviations get a half-student t prior with df = 3, mean = 0, and scale = 2.5
# these choices seem reasonable as default, so keep (for now)

# lasso prior
prior.lasso <- set_prior("lasso(df = 1, scale = 1)", class = "b") # defaults in brms
fit.lasso <- brm(yi|se(SEi) ~ . -vi - study + (1|study), data = scaled.training, prior = prior.lasso)

# horseshoe prior
prior.hs <- c(set_prior("horseshoe(df = 1, scale_global = 1, df_global = 1, scale_slab = 2, df_slab = 4)", class = "b")) # defaults in brms
fit.hs <- brm(yi|se(SEi) ~ . -vi - study + (1|study), data = scaled.training, prior = prior.hs)

##### Extract the relevant moderators -----
# Projection predictive variable selection does not seem to work because only 1 ES per study
fit.lasso.vs <- varsel(fit.lasso)

# Use CIs instead (for now)
# this is not ideal because 1) does not take into account correlations between moderators, and 2) which % CI to use is arbitrary
CI95 <- posterior_interval(fit.lasso, pars = paste0("X", 1:5), prob = 0.95)
sel.mod <- apply(CI95, 1, function(x){ !(sum(sign(x)) == 0)})
names(sel.mod) <- paste0("X", 1:5)
sel.mod
hpdi <- bayestestR::hdi(fit.hs)

# other output to return
summary(fit.lasso)

# density plots can already give a clear idea of which moderators are relevant
mcmc_areas(fit.lasso, pars = paste0("b_X", 1:5))
mcmc_areas(fit.hs, pars = paste0("b_X", 1:5))

##### Function to estimate the model -----
# data should be provided as df with yi, vi, SEi, and study variables (named as such)
brma <- function(data, moderators, prior = c("lasso", "horseshoe"), scale = TRUE, df_lasso = 1, scale_lasso = 1,
                 df_hs = 1, scale_global = 1, df_global = 1, scale_slab = 1, df_slab = 4,
                 iter = 2000, chains = 4){

  if(scale == TRUE){
    cols <- sapply(moderators, function(x) grep(x, colnames(data)))
    scaled.mods <- scale(data[, cols])
    data <- cbind.data.frame(data[, -cols], scaled.mods)
  }

  if(prior == "lasso"){
    prior <- set_prior(paste0("lasso(df = ", df_lasso, ", scale = ", scale_lasso, ")"),
                       class = "b")
  } else{
    prior <- set_prior(paste0("horseshoe(df = ", df_hs, ", scale_global = ", scale_global,
                              ", df_global = ", df_global, ", scale_slab = ", scale_slab,
                              ", df_slab = ", df_slab), class = "b")
  }

  formula <- paste("yi|se(SEi) ~ ", paste(moderators, collapse = "+"), " + (1|study)", sep = "")
  fit <- brm(formula, data = data, prior = prior, iter = iter, chains = chains)

  #TODO: possibly customize convergence warnings with tips on improving convergence given the possibilities within the brma function

  return(fit)

}

fit <- brma(simdat$training, moderators = paste0("X", 1:5), prior = "lasso")

##### Post processing functions to include -----
## density plots can already give a clear idea of which moderators are relevant
mcmc_areas(fit, pars = paste0("b_X", 1:5))

## variable selection based on CI
select_CI <- function(fit, moderators, prob){
  CI <- posterior_interval(fit, pars = moderators, prob = prob)

  sel <- apply(CI, 1, function(x){
    ifelse(x[1] < 0 & x[2] > 0, "no", "yes")
  })

  names(sel) <- moderators

  return(sel)
}

select_CI(fit = fit, moderators = paste0("X", 1:5), prob = 0.80)

## summary
summary(fit)






