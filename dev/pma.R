library(pema)
library(metafor)
library(glmnet)
ntau <- 10


df <- pema::simulate_smd(k_train = 1000, es = 2)

df_train <- df$training

pma <- function(x, y, vi, ntau = 10, ...){
  mf_nomod <- rma(y, vi)
  tau2_max <- mf_nomod$tau2
  mf_allmod <- rma(yi = y, vi = vi, mods = x)
  tau2_min <- mf_allmod$tau2
  tau_seq <- seq(from = tau2_min, to = tau2_max, length.out = ntau)
  wts <- 1/(df_train$vi + tau_seq[1])
  res <- vector(mode = "list", length = ntau)
  res[[1]] <- cv.glmnet(x = X, y = y, weights = wts, keep = TRUE)
  for(i in 2:ntau){
    wts <- 1/(df_train$vi + tau_seq[i])
    res[[i]] <- cv.glmnet(x = X, y = y, weights = wts, foldid = res[[1]]$foldid)
  }
  err_min <- sapply(res, function(i){min(i$cvm)})
  tau2_cv <- tau_seq[which.min(err_min)]
  list(
    rma = mf_nomod,
    lasso = res[[which.min(err_min)]],
    tau2_cv = tau2_cv
  )
}

m1 <- pma(x = as.matrix(df_train[, -c(1:2)]), y = df_train$yi, vi = df_train$vi)
m1$rma
m1$lasso
coef(m1$lasso)
m3 <- rma(yi = df_train$yi, vi = df_train$vi, mods = df_train[,-c(1:2)])
m3
mf_nomod <- rma(df_train$yi, df_train$vi)
tau2_max <- mf_nomod$tau2

mf_allmod <- rma(df_train$yi, df_train$vi, mods = df_train[, -c(1:2)])
tau2_min <- mf_allmod$tau2

tau_seq <- seq(from = tau2_min, to = tau2_max, length.out = ntau)

y <- df_train$yi
X <- as.matrix(df_train[, -c(1:2)])
wts <- 1/(df_train$vi + tau_seq[1])
res1 <- cv.glmnet(x = X, y = y, weights = wts, keep = TRUE)
res <- lapply(tau_seq, function(t2){
  wts <- 1/(df_train$vi + t2)
  cv.glmnet(x = X, y = y, weights = wts, foldid = res1$foldid)
})
errs <- sapply(seq_along(res), function(i){c(which.min(res[[i]]$cvm), min(res[[i]]$cvm))})
tmp <- which.min(errs[2,])
c(tau = tmp, lambda = errs[1, tmp])
})
tmp
sapply(seq_along(res), function(i){res[[i]]$cvm[tmp[2,i]]})
sapply(tmp[2,], function(i){c(which.min(res[[i]]$cvm), min(res[[i]]$cvm))})




X <- as.matrix(iris[,1:3])
y <- iris$Petal.Width
res <- replicate(100, {
  res1 <- cv.glmnet(x = X, y = y)
  res2 <- cv.glmnet(x = X, y = y, weights = rep(1, length(y)))
  c(which.min(res1$cvm), which.min(res2$cvm))
})
