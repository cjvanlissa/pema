library(metafor)
library(metaSEM)

::meta()

### calculate log risk ratios and corresponding sampling variances
dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

### fit a random-effects model using the log risk ratios and variances as input
### note: method="REML" is the default, so one could leave this out
rma(yi, vi, data=dat, method="REML")

tmp <- meta(y = dat[, "yi", drop = FALSE], v = dat[, "vi", drop = FALSE])
summary(mxRun(tmp$mx.model))
tmp$mx.model$Inter$values

X <- model.matrix(~. , dat[c("year", "alloc")])[, -1]

tmp <- meta(y = dat[, "yi", drop = FALSE], v = dat[, "vi", drop = FALSE], x = X, run= TRUE)
summary(mxTryHard(mxRun(tmp$mx.model)))
tmp$mx.model$Inter$values

tmp$

summary(tmp)

# install.packages("devtools")
devtools::install_github("trbrick/mxregsem")

install.packages("regsem")
library(regsem)
library(mxregsem)
?mxRegularizeLASSO
mxregsem::mxRegularizeLASSO()
tmp <- OpenMx::mxModel(tmp$mx.model, mxregsem::mxRegularizeLASSO(what = paste0("Slope1_", 1:3), name = NULL))
fit1 <- regsem(tmp$mx.model, lambda=0.05, type="lasso")
