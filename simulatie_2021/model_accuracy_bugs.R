library(pema)
library(brms)

#load in list of brma_models
#you can optionally run the models yourself (line 70 - 101)
load('./simulatie_2021/brma_models.RData')

simdat <- simulate_smd()
fit.lasso <- brma(yi~., data = simdat$training, method = "lasso")

source('./simulatie_2021/model_accuracy.R')

# Preparation of pred_brma to be run without model_accuracy dependency
pred_brma <- function(x, data = NULL, ...){
  # Prepare data ------------------------------------------------------------
  if(is.null(data)){
    fit <- x #just added this line to make sure it runs without model_accuracy() dependency
    X <- fit$X
  } else {
    df <- data
    if(!is.null(x[["vi_column"]])){
      df[[x$vi_column]] <- NULL
    }
    if(!is.null(x[["study_column"]])){
      df[[x[["study_column"]]]] <- NULL
    }
    mf <- call("model.frame",
               formula = x$formula,
               data = df,
               drop.unused.levels = TRUE)
    mf <- eval(mf, parent.frame())
    X <- mf[, -1, drop = FALSE]
  }

  # Prepare coefs -----------------------------------------------------------

  coefs <- rstan::summary(x$fit)$summary[, "mean"]
  int <- coefs[names(coefs) == "Intercept"]
  coefs <- coefs[startsWith(names(coefs), "betas[")]

  # Produce prediction ------------------------------------------------------
  int + rowSums(X * outer(rep.int(1L, nrow(X)), coefs))
}

#this works fine
out.pred1 <- pred_brma(fit.lasso)

#also when applied to list
out.pred.list <- lapply(brma_models, pred_brma)

#model_accuracy() however returns error: 'Error in df[[x$vi_column]] <- NULL : replacement has length zero'
model_accuracy(fit = fit.lasso,
  newdata = as.matrix(simdat$training[, -c(1:2)]),
  observed = simdat[[1]]$yi
)

#this error message seems to come from pred_brma. This is weird cause it works fine when ran independent of
#model_accuracy(). The error comes from this part:
if(!is.null(x[["vi_column"]])){
  df[[x$vi_column]] <- NULL
}

#the code gives same error if the code is changed to
if(!is.null(x[["vi_column"]])){
  df[[x[["vi_column"]]]] <- NULL
}



#OPTIONALLY run BRMA models yourself
library(parallel)
library(pema)
library(brms)

#make clusters
cores <- as.integer(floor(detectCores() / 2))
cl <- makeCluster(cores)

#export to clusters
clusterEvalQ(cl, library(pema))
clusterEvalQ(cl, library(brms))

#create list (length2) of simulated datasets
simulated_data <- list(pema::simulate_smd(), pema::simulate_smd())

#brma for simulation
brma_for_sim <- function(...){
  args <- as.list(match.call()[-1])
  res <- try(do.call(brma, args), silent = TRUE)
}

#export to clusters
clusterExport(cl, "brma_for_sim")

#run models
brma_models <- parLapply(cl, simulated_data, function(data) {
  args <- list(formula = force(as.formula(paste0("yi~", paste0(names(data$training)[-c(1:2)], collapse = "+")))),
               data = data$training,
               method = 'lasso')
  do.call(brma_for_sim, args)
})
