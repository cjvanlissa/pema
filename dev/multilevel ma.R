# A Comparison of the rma() and the lm(), lme(), and lmer() Functions
# In essence, the commonly used meta-analytic models are just special cases of
# the linear (mixed-effects) model with the only peculiar aspect that the
# variances of the error terms (i.e., the sampling variances) are known. Not
# surprisingly, the question therefore comes up occasionally why one cannot use
# the lm(), lme(), and lmer() functions for conducting meta-analyses. After
# all, the functions can be used to fit linear (mixed-effects) models and both
# functions allow the user to specify the sampling variances via the weights
# argument. So it seems that one should also be able to fit meta-analytic
# models with these functions. Below, I describe and illustrate how the models
# fitted via the lm(), lme(), and lmer() functions differ from the models
# fitted by the rma() function and why the those functions are therefore not
# suitable for fitting meta-analytic models.

# Fixed-Effects Model
# Let's start with a fixed-effects model. As an example, consider the data from
# the meta-analysis by Molloy et al. (2014) on the relationship between
# conscientiousness and medication adherence. For each study, we can compute
# the r-to-z transformed correlation coefficient and corresponding sampling
# variance with:

library(metafor)
dat <- escalc(measure="ZCOR", ri=ri, ni=ni, data=dat.molloy2014)
yi <- dat$yi
vi <- dat$vi
X <- dat[5:10]

# The yi values are the r-to-z transformed correlations and the vi values the
# corresponding sampling variances.

# We can now fit a fixed-effects model to these data with:

res.fe <- rma(yi, vi, data=dat, method="FE")
res.fe

# Now let's try to fit the same model with the lm() function by specifying the
# inverse of the sampling variances as weights:

res.lm <- lm(yi ~ 1, weights = 1/vi, data=dat)
summary(res.lm)

# Two things are of note here:

# The estimated intercept (0.12518) is exactly the same as the model
# coefficient obtained earlier with the rma() function (the value reported by
# the rma() function is rounded to 4 decimals, but that can be changed with
# print(res.fe, digits=5) yielding the same value).
# The standard error (of the estimated intercept) is different (and hence, the
# test statistic and p-value also differs).
# The reason for this discrepancy is that the model fitted by the lm() function
# assumes that the weights are not exactly known, but only up to a
# proportionality constant (namely s2e, that is, the error or residual
# variance), which is estimated and then factored into the weights.

# This can be demonstrated by extracting the estimated error variance from the
# lm object, multiplying the sampling variances by that value, and re-fitting
# the model with the rma() function:

rma(yi, vi*summary(res.lm)$sigma^2, data=dat, method="FE")

# Now these results are exactly the same as those obtained by the lm()
# function. Note that multiplying the sampling variances by some constant value
# does not affect the model estimate, but only the corresponding standard
# error.

# Therefore, if we want to obtain the same standard error from lm() as from
# rma(), we must adjust the standard error by ^se
# (i.e., the square-root of the estimated error variance or what is called the
# "residual standard error" in the output from the lm() function):

coef(summary(res.fe))

coef(summary(res.lm))[1,2] / summary(res.lm)$sigma

# Random-Effects Model
# A random-effects model can be fitted to the same data using the rma() function with:

res.re <- rma(yi, vi, data=dat)
res.re

# Now let's try to fit the same model with the lme() function from the nlme
# package. Here, the weights argument is used to specify a variance function
# with fixed variances (note: the name of the weights argument is a bit of a
# misnomer for the lme() function, as one does not use it to specify the
# weights – as in the lm() function – but a variance function structure). In
# addition, we need to add the random effects for each study via the random
# argument.

library(nlme)
dat$study <- 1:nrow(dat)
res.lme <- lme(yi ~ 1, random = ~ 1 | study, weights = varFixed(~ vi), data=dat)
summary(res.lme)

# Alternatively, we could use the lmer() function from the lme4 package. The syntax would be:

library(lme4)
res.lmer <- lmer(yi ~ 1 + (1 | study), weights = 1/vi, data=dat,
                 control=lmerControl(check.nobs.vs.nlev="ignore", check.nobs.vs.nRE="ignore"))
summary(res.lmer)

# Note that some of the default checks need to be switched off before lmer()
# will go through with fitting this model.

# A couple things are of note here:

# The estimated intercept differs from the model coefficient obtained earlier
# with the rma() function.
# The standard error is also different (and hence, the test statistic and
# p-value also differs).
# The estimated standard deviation of the random effects also differs (0.0682
# for lme() and lmer() compared to 0.0901 for rma()).
# These discrepancies are due to the exact same reason described earlier. The
# lme() and lmer() functions assume that the sampling variances are not exactly
# known, but again just up to a proportionality constant, namely the residual
# variance.

# To illustrate this, we can again factor in that constant into the sampling
# variances and refit the model with rma():

rma(yi, vi*res.lme$sigma^2, data=dat)

# Now these results match exactly what was obtained with the lme() function.

# In terms of the underlying models, we can also easily describe the difference
# as follows. The underlying model fitted by the rma() function is:
# yi = µ + ui + ei,
# where ui ~ N(0,t2) and ei ~ N(0,vi), where the vi
# values are the known sampling variances. On the other hand, the model fitted
# by the lme() and lmer() functions is:
# yi = µ + ui + ei,
# where ui ~ N(0, t2) and ei ~ N(0, s2e*vi), where
# s2e is that proportionality constant.

# Meta-Regression
# The examples above show models without any covariates/predictors/moderators
# (i.e., the models only include an intercept term). However, the exact same
# discrepancy will be found when including covariates into the models. In this
# context, Thompson and Sharp (1999) also describe the difference discussed
# above. They denote the model fitted by lm() as a model that accounts for
# potential heterogeneity by a multiplicative factor (i.e., the residual
# variance), while the random-effects model fitted via the rma() function is a
# model that accounts for potential heterogeneity by an additive factor. They
# note that the "rationale for using a multiplicative factor for variance       !
# inflation is weak" (p. 2705) and clearly recommend to use a model with an     !
# additive heterogeneity component.

# The possibility of using both an additive and multiplicative component
# simultaneously was not discussed by Thompson and Sharp (1999). This is in
# fact what the model fitted with the lme() function does. However, as far as I
# am aware of, no motivation supporting such a model has ever been given. In    !
# fact, it is my impression that those who use the lme() and lmer() functions   !
# for meta-analytic purposes are typically not aware of the issues discussed    !
# above.

# More Complex Models
# Note that the same issues apply when fitting more complex meta-analytic
# models (e.g., multivariate/multilevel models). The rma.mv() function in the
# metafor package can be used to fit such models (see the analysis examples
# section for some illustrations). Using the lme() and lmer() functions to fit
# such models will again lead to a mix of additive and multiplicative variance
# components, which is not how multivariate/multilevel meta-analytic models are
# typically defined.

# S-Plus Version of lme()
# The nlme package was developed by José Pinheiro and Douglas Bates for both R
# and S-Plus. Interestingly, the S-Plus version has a special control argument  !
# that allows the user to set the error variance to a fixed constant. By        !
# setting this to 1, one can fit the exact same model as the rma() does:

res.lme <- lme(yi ~ 1, random = ~ 1 | study, weights = varFixed(~ vi), control=lmeControl(sigma = 1), data=dat)
summary(res.lme)
res.re

# These are the exact same results as obtained earlier with the rma() function.
# Unfortunately, the R version of the nlme package does not provide this
# functionality.

# Update: The R version of the nlme package does allow the use of the
# lmeControl(sigma = 1) control argument (this was added in version 3.1-123,
# which was released 2016-01-17). However, using this does not yield the same
# results as obtained above (the results are close but not the same). The
# results given above have also been compared to those obtained with Stata
# (metareg command) and SAS (proc mixed) and are definitely correct. This issue
# has also been raised on the R-sig-ME mailing list (see here). See also this
# blog post by James Pustejovsky for some further comparisons between lme() and
# the rma() and rma.mv() functions.

# Summary
# The models fitted by the rma() function assume that the sampling variances
# are known. The models fitted by the lm(), lme(), and lmer() functions assume
# that the sampling variances are known only up to a proportionality constant.
# These are different models than typically used in meta-analyses.

# References
# Molloy, G. J., O'Carroll, R. E., & Ferguson, E. (2014). Conscientiousness and
# medication adherence: A meta-analysis. Annals of Behavioral Medicine, 47(1),
# 92–101.

# Thompson, S. G., & Sharp, S. J. (1999). Explaining heterogeneity in
# meta-analysis: A comparison of methods. Statistics in Medicine, 18(20),
# 2693–2708.

# tips/rma_vs_lm_lme_lmer.txt · Last modified: 2021/04/16 07:50 by Wolfgang
# Viechtbauer
# Page Tools
# Except where otherwise noted, content on this wiki is licensed under the
# following license: CC Attribution-Noncommercial-Share Alike 4.0 International


# https://www.jepusto.com/bug-in-nlme-with-fixed-sigma/



# Now with moderators -----------------------------------------------------


res.lme <- lme(as.formula(paste0("yi ~ 1 + ", paste0(names(X), collapse = " + "))), random = ~ 1 | study, weights = varFixed(~ vi), control=lmeControl(sigma = 1), data=dat)
summary(res.lme)
res.re


library(glmmLasso)
glmmLasso(fix=yi ~ 1, rnd = ~ 1 | study, data = dat, lambda)
dat$studfac <- as.factor(dat$study)
res.glml <- glmmLasso(fix=yi ~ 1 + as.factor(controls) + as.factor(design) + as.factor(a_measure)+ as.factor(c_measure) + meanage + quality, rnd = list(studfac = ~ 1), data = dat, lambda = 1)

round(cbind(res.glml$coefficients, res.lme$coefficients$fixed), 3)

dat$int <- 1L
res.glml <- glmmLasso(fix=yi ~ -1 + int, rnd = list(studfac = ~ 1), data = dat, lambda = 1)
