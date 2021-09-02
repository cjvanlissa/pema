# # Install packages; skip if already installed -----------------------------
# install.packages('devtools')                #to download pema from source
# devtools::install_github("cjvanlissa/pema") #for brma and simulate_smd functions
# install.packages("metaforest")              #for random effects metaforest and rma function
# install.packages("parallel")                #for parallel computing
# install.packages("sn")                      #for rsn function which allows skewed distributions
# install.packages("brms")                    #for package Rcpp that interfaces C++ code into R
# install.packages("rstan")                   #to compile C++ in Rstudio
# install.packages("rstantools")              #actually not sure if necessary
# install.packages("bayestestR")              #actually not sure if necessary

# Load packages
packagenames <- c('pema', 'metaforest', 'parallel', 'sn', 'brms', 'rstan', 'bayestestR', 'foreach')
lapply(packagenames, library, character.only = TRUE)

# Check package versions
# Which packages require correct versions and which versions are they?--------------------------------------------------
versions <- c(
  compareVersion(as.character(packageVersion("metaforest")), "0.1.3"),
  compareVersion(as.character(packageVersion("pema")), "0.1.0.9000"),
  compareVersion(as.character(packageVersion("rstan")), "2.21.2"))
if(!all(versions == 0)) stop("Using the incorrect version of one or more packages.")

# Load simulation functions from source -----------------------------------
source("./simulatie_2021/functions for simulation/simulate_smd.R") #this function is necessary to load as it seems to be
#different than the function simulate_smd in pema.
source("./simulatie_2021/functions for simulation/model_accuracy.R")
# source("./simulatie_2021/functions for simulation/mf_sim.R")
# source("./simulatie_2021/functions for simulation/rma_sim.R")
# source("./simulatie_2021/functions for simulation/rma_for_sim.R")
# source("./simulatie_2021/functions for simulation/brma_for_sim.R") #contains both brma_for_sim and brma_sim functions

#set conditions for simulation
hyper_parameters<-list(
  #Number of datasets per condition
  ndataset=1:100,
  #Number of studies per dataset, normally distributed with mean n and sd n/3
  k_train=c(20, 40, 100),
  #Average n per study (k)
  mean_n=c(40, 80, 160),
  #Effect size
  es = c(0, .2, .5, .8),
  #Residual heterogeneity
  tau2=c(0, .04, .1),
  # Slant parameter alpha
  alpha_tau = c(0),
  alpha_mod = c(0, 2, 10),
  #Study-level moderators
  moderators= c(1, 2, 5),
  # This model includes effects of the first 10 moderators. The next 10 moderators are irrelevant.
  model = c("es * x[, 1]",
            "es * x[, 1] + es * x[, 2] + es * (x[, 1] * x[, 2])",
            "es * x[, 1] + es * (x[, 1] ^ 2) + es * (x[, 1] ^ 3)",
            "es * exp(x[, 1])"),
  distribution = "sn"
)

# Create grid with simulation parameters and save it
summarydata <- expand.grid(hyper_parameters, stringsAsFactors = FALSE)

# Add moderators so there's always sufficient predictors
summarydata$moderators <-  summarydata$moderators + c(1, 2, 1, 1)[factor(summarydata$model, levels = c("es * x[, 1]",
                                                                                                       "es * x[, 1] + es * x[, 2] + es * (x[, 1] * x[, 2])",
                                                                                                     "es * x[, 1] + es * (x[, 1] ^ 2) + es * (x[, 1] ^ 3)",
                                                                                                       "es * exp(x[, 1])"))]
saveRDS(summarydata, file = "summarydata.RData")
#summarydata<-readRDS("summarydata_2017-09-22_16_06.RData")


#Start clusters, load relevant packages and export functions to clusters
# no_cores <- 12#detectCores()
# cl <- makeCluster(no_cores)

library(doSNOW)

# clusterEvalQ(cl, library(metaforest))
# clusterEvalQ(cl, library(data.table))
# clusterEvalQ(cl, library(sn))
# clusterEvalQ(cl, library(pema))
# clusterEvalQ(cl, library(brms))
# clusterEvalQ(cl, library(rstan))
# clusterEvalQ(cl, library(bayestestR))
# clusterExport(cl, "simulate_smd")
# clusterExport(cl, "model_accuracy")

set.seed(78326)
library(foreach)
cl<-makeCluster(40) #change the 2 to your number of CPU cores
registerDoSNOW(cl)

# f <- list.files("./simulatie_2021/results", full.names = TRUE)
# file.remove(f)
out_sel <- rep(NA, max(summarydata$moderators))
foreach(rownum = 1:nrow(summarydata), .packages = c("pema", "sn", "metaforest", "rstan", "bayestestR", "brms")) %dopar% {
  data <- do.call(simulate_smd, as.list(summarydata[rownum, -1, drop = FALSE]))
  models <- tryCatch(MetaForest(yi ~., data = data$training, num.trees = 100, whichweights = "random"),
                     error = function(e){ NULL })
  if(!is.null(models)){
    fits <- c(
      rownum,
      model_accuracy(
        fit = models,
        newdata = subset(data[[1]], select = -c(1, 2)),
        observed = data[[1]]$yi
      ),
      model_accuracy(
        fit = models,
        newdata = subset(data[[2]], select = -1),
        observed = data[[2]]$yi,
        ymean = mean(data[[1]]$yi)
      ),
      tau2 = models$rma_after$tau2
    )
    #importance_pvalues(models$forest, method = "janitza")[, 2] < .05
    out <- c(rownum, out_sel)
    out[2:(length(models$forest$variable.importance)+1)] <- models$forest$variable.importance > 0
  } else {
    fits <- c(
      rownum,
      rep(NA, 7))
    out <- c(rownum, out_sel)
  }
  write.table(x = t(fits), file = sprintf("./simulatie_2021/results/mf_%d.txt" , Sys.getpid()), sep = "\t", append = TRUE, row.names = FALSE, col.names = FALSE)
  write.table(x = t(out), file = sprintf("./simulatie_2021/results/mf_sel_%d.txt" , Sys.getpid()), sep = "\t", append = TRUE, row.names = FALSE, col.names = FALSE)


# Conduct rma models ------------------------------------------------------

  models <- tryCatch(rma(yi = data$training$yi,
                         vi = data$training$vi,
                         mods = as.matrix(data$training[,-c(1:2)])),
                     error = function(e){ NULL })
  if(!is.null(models)){
    fits <-
      c(rownum,
        model_accuracy(
          fit = models,
          newdata = as.matrix(data$training[, -c(1:2)]),
          observed = data[[1]]$yi
        ),
        model_accuracy(
          fit = models,
          newdata = as.matrix(data$testing[, -1]),
          observed = data$testing$yi,
          ymean = mean(data$training$yi)
        ),
        tau2 = models$tau2,
        beta1 = models$b[2,1]
      )
    out <- c(rownum, out_sel)
    out[2:length(models$pval)] <- models$pval[-1] < .05
  } else {
    fits <- c(
      rownum,
      rep(NA, 8))
    out <- c(rownum, out_sel)
  }
  write.table(x = t(fits), file = sprintf("./simulatie_2021/results/rma_%d.txt" , Sys.getpid()), sep = "\t", append = TRUE, row.names = FALSE, col.names = FALSE)
  write.table(x = t(out), file = sprintf("./simulatie_2021/results/rma_sel_%d.txt" , Sys.getpid()), sep = "\t", append = TRUE, row.names = FALSE, col.names = FALSE)




# Conduct pema models ------------------------------------------------------

  models <- tryCatch(brma(formula = yi~.,
                          data = data$training,
                          method = "lasso",
                          chains = 1,
                          iter = 6000,
                          refresh = -1,
                          verbose = FALSE),
                     error = function(e){ NULL })
  if(!is.null(models)){
    summ <- summary(models$fit)$summary
    fits <-
      c(rownum,
        model_accuracy(
          fit = models,
          newdata = data$training[, -1],
          observed = data[[1]]$yi
        ),
        model_accuracy(
          fit = models,
          newdata = data$testing,
          observed = data$testing$yi,
          ymean = mean(data$training$yi)
        ),
        tau2 = mean(models$fit@sim$samples[[1]]$`sd_1[1]`),
        beta1 = mean(models$fit@sim$samples[[1]]$`betas[1]`),
        divergent = {
          np = nuts_params(models$fit); sum(np$Parameter == "divergent__" & np$Value == 1)
        },
        maxRhat = max(summ[, "Rhat"], na.rm = TRUE),
        minNeff = min(summ[, "n_eff"], na.rm = TRUE)
      )
    outCI <- outHDI <- out_sel
    CI95 <- summary(models$fit)$summary[paste0("betas[", 1:(ncol(models$X)-1), "]"), c("2.5%", "97.5%"), drop = FALSE]
    selected_CI95 <- apply(CI95, 1, function(x){ !(sum(sign(x)) == 0)})
    outCI[1:length(selected_CI95)] <- selected_CI95
    hpdi <- bayestestR::hdi(models$fit, parameters = "betas")
    hpdi <- cbind(hpdi$CI_low, hpdi$CI_high)
    outHDI[1:length(selected_CI95)] <- apply(CI95, 1, function(x){ !(sum(sign(x)) == 0)})
    out <- c(rownum, outCI, outHDI)
  } else {
    fits <- c(
      rownum,
      rep(NA, 12))
    out <- c(rownum, out_sel, out_sel)
  }
  write.table(x = t(fits), file = sprintf("./simulatie_2021/results/las_%d.txt" , Sys.getpid()), sep = "\t", append = TRUE, row.names = FALSE, col.names = FALSE)
  write.table(x = t(out), file = sprintf("./simulatie_2021/results/las_sel_%d.txt" , Sys.getpid()), sep = "\t", append = TRUE, row.names = FALSE, col.names = FALSE)

  models <- tryCatch(brma(formula = yi~.,
                          data = data$training,
                          method = "hs",
                          chains = 1,
                          iter = 6000,
                          refresh = -1,
                          verbose = FALSE),
                     error = function(e){ NULL })
  if(!is.null(models)){
    summ <- summary(models$fit)$summary
    fits <-
      c(rownum,
        model_accuracy(
          fit = models,
          newdata = data$training[, -1],
          observed = data[[1]]$yi
        ),
        model_accuracy(
          fit = models,
          newdata = data$testing,
          observed = data$testing$yi,
          ymean = mean(data$training$yi)
        ),
        tau2 = mean(models$fit@sim$samples[[1]]$`sd_1[1]`),
        beta1 = mean(models$fit@sim$samples[[1]]$`betas[1]`),
        divergent = {
          np = nuts_params(models$fit); sum(np$Parameter == "divergent__" & np$Value == 1)
        },
        maxRhat = max(summ[, "Rhat"], na.rm = TRUE),
        minNeff = min(summ[, "n_eff"], na.rm = TRUE)
      )
    outCI <- outHDI <- out_sel
    CI95 <- summary(models$fit)$summary[paste0("betas[", 1:(ncol(models$X)-1), "]"), c("2.5%", "97.5%"), drop = FALSE]
    selected_CI95 <- apply(CI95, 1, function(x){ !(sum(sign(x)) == 0)})
    outCI[1:length(selected_CI95)] <- selected_CI95
    hpdi <- bayestestR::hdi(models$fit, parameters = "betas")
    hpdi <- cbind(hpdi$CI_low, hpdi$CI_high)
    outHDI[1:length(selected_CI95)] <- apply(CI95, 1, function(x){ !(sum(sign(x)) == 0)})
    out <- c(rownum, outCI, outHDI)
  } else {
    fits <- c(
      rownum,
      rep(NA, 12))
    out <- c(rownum, out_sel, out_sel)
  }
  write.table(x = t(fits), file = sprintf("./simulatie_2021/results/hs_%d.txt" , Sys.getpid()), sep = "\t", append = TRUE, row.names = FALSE, col.names = FALSE)
  write.table(x = t(out), file = sprintf("./simulatie_2021/results/hs_sel_%d.txt" , Sys.getpid()), sep = "\t", append = TRUE, row.names = FALSE, col.names = FALSE)
}

#Close cluster
stopCluster(cl)


# End of simulation -------------------------------------------------------

stop("End of simulation")

# Merge files -------------------------------------------------------------
library(data.table)
res <- as.data.table(readRDS(list.files(pattern = "summarydata")))

model_names <- unique(gsub("_fits_\\d+.RData", "", list.files(pattern="fits", all.files=FALSE, full.names=FALSE)))

# fit_vars <- c("_train_r2", "_train_mse", "_train_r", "_test_r2", "_test_mse", "_test_r", "_tau2")
# fit_vars <- paste0(unlist(lapply(model_names, rep, length(fit_vars))), fit_vars)
fit_vars <-
list(
  hs = c(
    "hs_train_r2",
    "hs_train_mse",
    "hs_train_r",
    "hs_test_r2",
    "hs_test_mse",
    "hs_test_r",
    "hs_tau2",
    "hs_beta1",
    "hs_div",
    "hs_rhat",
    "hs_neff"
  ),
  las = c(
    "lasso_train_r2",
    "lasso_train_mse",
    "lasso_train_r",
    "lasso_test_r2",
    "lasso_test_mse",
    "lasso_test_r",
    "lasso_tau2",
    "lasso_beta1",
    "lasso_div",
    "lasso_rhat",
    "lasso_neff"
  ),
  mf = c(
    "mf_r_train_r2",
    "mf_r_train_mse",
    "mf_r_train_r",
    "mf_r_test_r2",
    "mf_r_test_mse",
    "mf_r_test_r",
    "mf_r_tau2"
  ),
  rma = c(
    "rma_train_r2",
    "rma_train_mse",
    "rma_train_r",
    "rma_test_r2",
    "rma_test_mse",
    "rma_test_r",
    "rma_tau2",
    "rma_beta1"
  )
)
fit_vars <- lapply(fit_vars, function(i){c("rownum", i)})

for(mod in names(fit_vars)){
  f <- list.files("./simulatie_2021/results", pattern = paste0("^", mod, "_\\d"), full.names = TRUE)
  tab <- fread_many(f)
  setorderv(tab, cols = "V1", order=1L, na.last=FALSE)
  if(!(tab$V1[1] == 1 & tail(tab$V1, 1) == nrow(res) & length(unique(tab$V1)) == nrow(res))){
    stop()
  }
  names(tab) <- fit_vars[[mod]]
  tab[, "rownum" := NULL]
  res <- cbind(res, tab)
}

for(mod in names(fit_vars)){
  f <- list.files("./simulatie_2021/results", pattern = paste0("^", mod, "_sel_\\d"), full.names = TRUE)
  tab <- fread_many(f)
  setorderv(tab, cols = "V1", order=1L, na.last=FALSE)
  if(!(tab$V1[1] == 1 & tail(tab$V1, 1) == nrow(res) & length(unique(tab$V1)) == nrow(res))){
    stop()
  }
  tab[, "V1" := NULL]
  if(mod %in% c("hs", "las")){
    names(tab) <- paste0(mod, rep(c("_ci_", "_hdi_"), each = ncol(tab)/2), "sel_", 1:(ncol(tab)/2))
  } else {
    names(tab) <- paste0(mod, "_sel_", 1:ncol(tab))
  }
  res <- cbind(res, tab)
}

fwrite(res, file.path("simulatie_2021", paste0("sim_results_", Sys.Date(), ".csv")))
saveRDS(res, file.path("simulatie_2021", paste0("sim_results_", Sys.Date(), ".RData")))
f <- list.files("./simulatie_2021/results", full.names = TRUE)
file.remove(f)

# Delete alle variabelen waar jij niets mee hoeft
res[, grep("(_train_r2|_r|_mse)$", names(res), value = TRUE) := NULL]


# Tot hier is alles gekuisd -----------------------------------------------



library(ggplot2)
library(data.table)

vars <- c("_train_mse", "_train_mse", "_train_r", "_test_mse", "_test_mse", "_test_r")






?paste()


###################################
# Predictive performance analyses #
###################################

analyzedat <- data[ , .SD, .SDcols=c("k_studies", "mean_study_n", "es", "residual_heterogeneity", "moderators", "model", paste0(newnames, "_train_r2"))]
analyzedat[, c(1:6):=lapply(.SD, factor), .SDcols=c(1:6)]

sink("output.txt")
tmp<-rbind(round(analyzedat[model != "1", lapply(.SD, mean), .SDcols=c(7:ncol(analyzedat))], 2),
           round(analyzedat[model !="1", lapply(.SD, sd), .SDcols=c(7:ncol(analyzedat))], 2))
tmp<-data.frame(names(tmp), t(tmp))
apply(tmp, 1, function(x){paste(x, collapse=", ")})
sink()

#Predictive performance on testing data
analyzedat <- data[ , .SD, .SDcols=c("k_studies", "mean_study_n", "es", "residual_heterogeneity", "moderators", "model", paste0(newnames, "_test_r2"))]
analyzedat[, c(1:6):=lapply(.SD, factor), .SDcols=c(1:6)]


analyzedat[ , paste(newnames[1], newnames[-1], sep = "_") :=
              lapply(paste0(newnames[-1], "_test_r2"), function(x){
                analyzedat[[paste0(newnames[1], "_test_r2")]] - analyzedat[[x]]
              })]

tmp<-rbind(round(analyzedat[model=="1", lapply(.SD, mean), .SDcols=c(7:ncol(analyzedat))], 3),
           round(analyzedat[model=="1", lapply(.SD, sd), .SDcols=c(7:ncol(analyzedat))], 3))
tmp<-data.frame(names(tmp), t(tmp))
sink("output.txt", append = TRUE)
apply(tmp, 1, function(x){paste(x, collapse=", ")})
sink()

tmp<-rbind(round(analyzedat[model!="1", lapply(.SD, mean), .SDcols=c(7:ncol(analyzedat))], 3),
           round(analyzedat[model!="1", lapply(.SD, sd), .SDcols=c(7:ncol(analyzedat))], 3))
tmp<-data.frame(names(tmp), t(tmp))
sink("output.txt", append = TRUE)
apply(tmp, 1, function(x){paste(x, collapse=", ")})
sink()

#Check for significant moderators
EtaSq<-function (x)
{
  anovaResults <- summary.aov(x)[[1]]
  anovaResultsNames <- rownames(anovaResults)
  SS <- anovaResults[,2] #SS effects and residuals
  k <- length(SS) - 1  # Number of factors
  ssResid <- SS[k + 1]  # Sum of Squares Residual
  ssTot <- sum(SS)  # Sum of Squares Total
  SS <- SS[1:k] # takes only the effect SS
  anovaResultsNames <- anovaResultsNames[1:k]
  etaSquared <- SS/ssTot # Should be the same as R^2 values
  partialEtaSquared <- SS/(SS + ssResid)
  res <- cbind(etaSquared, partialEtaSquared)
  colnames(res) <- c("Eta^2", "Partial Eta^2")
  rownames(res) <- anovaResultsNames
  return(res)
}

yvars<-c(paste0(newnames, "_test_r2"), paste(newnames[1], newnames[-1], sep = "_"))


withoutmodel1<-analyzedat[model!="1",]
anovas<-lapply(yvars, function(yvar){
  form<-paste(yvar, "(es + residual_heterogeneity + k_studies + mean_study_n + moderators + model) ^ 2", sep="~")
  thisaov<-aov(as.formula(form), data=withoutmodel1)
  thisetasq<-EtaSq(thisaov)[ , 2]
  list(thisaov, thisetasq)
})

etasqs<-data.frame(sapply(anovas, function(x){
  tempnums<-x[[2]]
  formatC(tempnums, 2, format="f")
}))

coef.table<-etasqs
names(coef.table)<-yvars

table.names<-row.names(coef.table)
table.names<-gsub("_studies", "", table.names)
table.names<-gsub("moderators", "M", table.names)
table.names<-gsub("residual_heterogeneity", "$\\tau^2$", table.names)
table.names<-gsub("es", "$\\beta$", table.names)
table.names<-gsub("mean_study_n", "$\\mean{n}$", table.names)
table.names<-gsub("model", "Model ", table.names)

row.names(coef.table)<-table.names
names(coef.table)<-c("MetaForest", "MF_Knowntau", "MF_Residualtau" ,"Fixed-effects",  "Unweighted", "metaCART", "$\\Delta_{MF-MK}$", "$\\Delta_{MF-MR}$", "$\\Delta_{MF-FX}$", "$\\Delta_{MF-UN}$",   "$\\Delta_{MF-CA}$")

study1.coef.table.withoutmodel1 <- coef.table
sink("output.txt", append = TRUE)
study1.coef.table.withoutmodel1
sink()

study1.coef.table.withoutmodel1 <- study1.coef.table.withoutmodel1[, -c(2,3, 7, 8)]


onlymodel1<-analyzedat[model=="1",]
anovas<-lapply(yvars, function(yvar){
  form<-paste(yvar, "(es + residual_heterogeneity + k_studies + mean_study_n + moderators) ^ 2", sep="~")
  thisaov<-aov(as.formula(form), data=onlymodel1)
  thisetasq<-EtaSq(thisaov)[ , 2]
  list(thisaov, thisetasq)
})

etasqs<-data.frame(sapply(anovas, function(x){
  tempnums<-x[[2]]
  formatC(tempnums, 2, format="f")
}))

coef.table<-etasqs
names(coef.table)<-yvars

table.names<-row.names(coef.table)
table.names<-gsub("_studies", "", table.names)
table.names<-gsub("moderators", "M", table.names)
table.names<-gsub("residual_heterogeneity", "$\\tau^2$", table.names)
table.names<-gsub("es", "$\\beta$", table.names)
table.names<-gsub("mean_study_n", "$\\mean{n}$", table.names)
table.names<-gsub("model", "Model ", table.names)

row.names(coef.table)<-table.names

names(coef.table)<-c("MetaForest", "MF_Knowntau", "MF_Residualtau" ,"Fixed-effects",  "Unweighted", "metaCART", "$\\Delta_{MF-MK}$", "$\\Delta_{MF-MR}$", "$\\Delta_{MF-FX}$", "$\\Delta_{MF-UN}$",   "$\\Delta_{MF-CA}$")

study1.coef.table.onlymodel1 <- coef.table
sink("output.txt", append = TRUE)
coef.table
sink()
study1.coef.table.onlymodel1 <- study1.coef.table.onlymodel1[, -c(2,3, 7, 8)]

library(xtable)

print(xtable(study1.coef.table.onlymodel1), file="table_study1_model1_full.tex",sanitize.text.function=function(x){x})

print(xtable(study1.coef.table.withoutmodel1), file="table_study1_model2_5_full.tex",sanitize.text.function=function(x){x})


##Figure out which predictors are most important per variable
as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}
sink("output.txt", append = TRUE)
lapply(study1.coef.table.withoutmodel1, function(x){
  tmp <- as.numeric.factor(x)
  names(tmp)<- c("k", "n", "beta", "tau2", "m", "model", "k:n", "k:beta", "k:tau2", "k:m", "k:model", "n:beta", "n:tau2", "n:m", "n:model", "beta:tau2", "beta:m", "beta:model", "tau2:m", "tau2:model", "m:model")
  sort(tmp, decreasing = TRUE)
})
sink()

##################
# Power analysis #
##################
analyzedat <- data[ , .SD, .SDcols=c("k_studies", "mean_study_n", "es", "residual_heterogeneity", "moderators", "model", "mf_test_r2", "ca_test_r2")]
analyzedat[, c(1:6):=lapply(.SD, factor), .SDcols=c(1:6)]

measure.vars <- names(analyzedat)[-c(1:6)]
grouping.vars <- quote(list(k_studies,
                            mean_study_n,
                            es,
                            residual_heterogeneity,
                            moderators, model))

power2<-analyzedat[,lapply(.SD, function(x){ifelse(quantile(x, probs = .2, na.rm = TRUE)>0, 1, 0)}),by=eval(grouping.vars), .SDcols=measure.vars]

plotdat<-power2

plotdat[, power := "Neither"]
plotdat[mf_test_r2 == 1 & ca_test_r2 == 0, power := "Only MetaForest"]
plotdat[mf_test_r2 == 0 & ca_test_r2 == 1, power := "Only metaCART"]
plotdat[mf_test_r2 == 1 & ca_test_r2 == 1, power := "Both"]
plotdat[, power := factor(power)]
plotdat$power <- ordered(plotdat$power, levels = c("Neither", "Only MetaForest", "Only metaCART", "Both"))
table(plotdat$power)
names(plotdat)[c(1:6)]<-c("k", "n", "es", "tau2", "M", "Model")

categorical<-c("k", "n", "es", "tau2", "M", "Model")

plotdat[, (categorical) := lapply(.SD, factor), .SDcols=categorical]
plotdat$tau2 <- factor(plotdat$tau2, labels = c("tau^2:~.00", "tau^2:~.04", "tau^2:~.28"))
plotdat$es <- factor(plotdat$es, labels = c("beta:~0.2", "beta:~0.5", "beta:~0.8"))
levels(plotdat$Model)<-c("a", "b", "c", "d", "e")
table(plotdat$power)
thisplot <- ggplot(plotdat, aes(x=n, y=k, fill = power)) +
  geom_raster(hjust = 0, vjust = 0)+
  facet_grid(Model+M ~ es+tau2, labeller = labeller(es = label_parsed, tau2 = label_parsed, Model=label_both, M = label_both))+
  scale_fill_manual(values=c("white", "grey50", "grey25", "black"))+
  scale_x_discrete(expand=c(0,0))+
  scale_y_discrete(expand=c(0,0))+
  labs(x=expression(bar(n)), fill="Power > .80")
#theme(legend.title = element_blank()) +
#thisplot+annotate("text", label = "p[mf<mc]", parse=TRUE)

ggsave("power_plot.pdf", plot=thisplot, device="pdf", width=210, height=297, units="mm")

#Percentage table
power <- analyzedat[!(model==1),lapply(.SD, function(x){sum(x > 0)/length(x)}),by=eval(grouping.vars), .SDcols=measure.vars]
powertmp <- analyzedat[model==1,lapply(.SD, function(x){sum(x < 0)/length(x)}),by=eval(grouping.vars), .SDcols=measure.vars]
power <- rbindlist(list(powertmp, power))

power <- power
power <- melt(power, measure.vars = c("mf_test_r2", "ca_test_r2"))

tmp <- dcast(power, model+moderators+residual_heterogeneity~variable+ es + k_studies + mean_study_n)
write.csv(tmp, "study1 power.csv")

###################################################
#               R-squared plots                   #
###################################################
includevars <- c("mf", "fx", "un", "ca")
plotdat <- data[ , .SD, .SDcols=c("k_studies", "mean_study_n", "es", "residual_heterogeneity", "moderators", "model", paste0(includevars, "_test_r2"))]
plotdat[, c(1:6):=lapply(.SD, factor), .SDcols=c(1:6)]

id.vars=names(plotdat)[1:6]
plotdat2<-melt.data.table(plotdat, id.vars=id.vars)
names(plotdat2)<-c("k", "n", "es", "tau2", "Moderators", "Model", "Algorithm", "R2")
categorical<-c("es", "tau2", "Moderators", "Model", "Algorithm")
plotdat2[, (categorical) := lapply(.SD, factor), .SDcols=categorical]
levels(plotdat2$Algorithm)
plotdat2$Algorithm <- ordered(plotdat2$Algorithm, levels = c("mf_test_r2", "fx_test_r2", "un_test_r2", "ca_test_r2"))#, "mk_test_r2", "mr_test_r2"
levels(plotdat2$Algorithm)<-c("MetaForest", "MetaForest FX", "MetaForest UN", "metaCART")#, "MetaForest KT", "MetaForest RT"
plotdat2$tau2 <- factor(plotdat2$tau2, labels = c("tau^2:~.00", "tau^2:~.04", "tau^2:~.28"))
plotdat2$es <- factor(plotdat2$es, labels = c("beta:~0.2", "beta:~0.5", "beta:~0.8"))
plotdat2$n <- factor(plotdat2$n, labels = c("bar(n):~40", "bar(n):~80", "bar(n):~160"))
levels(plotdat2$Model)<-c("a", "b", "c", "d", "e")
#plotdat2 <- plotdat2[!(Algorithm %in% c("MetaForest FX", "MetaForest UN")),]

#Interaction between beta:Model
plotdat3<-plotdat2

plotdat3$es <- factor(plotdat3$es, labels = c("0.2", "0.5", "0.8"))
plots <- ggplot(plotdat3, aes(x = es, y = R2, group=Algorithm, linetype = Algorithm, shape = Algorithm)) +
  stat_summary(fun.y = mean,
               geom = "point") +
  stat_summary(fun.y = mean,
               geom = "path") +
  #scale_x_discrete(expand = c(0,0))+
  #scale_x_continuous(breaks=c(.2, .5, .8), limits=c(.2,.8))+
  facet_grid(. ~ Model, labeller = labeller(Model=label_both))+
  theme_bw()+ labs(y=expression(R^2), x=expression(beta))+
  geom_hline(yintercept = 0)
#plots
ggsave("int_beta_model_all_algos.pdf", plot=plots, device="pdf")

#Interaction between beta:M
plots <- ggplot(plotdat3, aes(x = es, y = R2, group=Algorithm, linetype = Algorithm, shape = Algorithm)) +
  stat_summary(fun.y = mean,
               geom = "point") +
  stat_summary(fun.y = mean,
               geom = "path") +
  #scale_x_continuous(breaks=c(.2, .5, .8), limits=c(.2,.8))+
  facet_grid(. ~ Moderators, labeller = labeller(Moderators = label_both))+
  theme_bw()+ labs(y=expression(R^2), x=expression(beta))+
  geom_hline(yintercept = 0)
plots
ggsave("int_beta_moderators_all_algos.pdf", plot=plots, device="pdf")

#Interaction between k:Model
plotdat3<-plotdat2
#plotdat3$k <- as.numeric.factor(plotdat3$k)
plots <- ggplot(plotdat3, aes(x = k, y = R2, group=Algorithm, linetype = Algorithm, shape = Algorithm)) +
  stat_summary(fun.y = mean,
               geom = "point") +
  stat_summary(fun.y = mean,
               geom = "path") +
  #scale_x_continuous(breaks=c(20, 40, 80, 120), limits=c(20, 120))+
  facet_grid(. ~ Model, labeller = labeller(Model=label_both))+
  theme_bw()+ labs(y=expression(R^2))+
  geom_hline(yintercept = 0)
plots
ggsave("int_k_model_all_algos.pdf", plot=plots, device="pdf")

#Interaction between k:beta
plots <- ggplot(plotdat3, aes(x = k, y = R2, group=Algorithm, linetype = Algorithm, shape = Algorithm)) +
  stat_summary(fun.y = mean,
               geom = "point") +
  stat_summary(fun.y = mean,
               geom = "path") +
  #scale_x_continuous(breaks=c(20, 40, 80, 120), limits=c(20, 120))+
  facet_grid(. ~ es, labeller = labeller(es = label_parsed))+
  theme_bw()+ labs(y=expression(R^2))+
  geom_hline(yintercept = 0)
#plots
ggsave("int_k_beta_all_algos.pdf", plot=plots, device="pdf")


#Marginal effect of k
plots <- ggplot(plotdat3, aes(x = k, y = R2, group=Algorithm, linetype = Algorithm, shape = Algorithm)) +
  stat_summary(fun.y = mean,
               geom = "point") +
  stat_summary(fun.y = mean,
               geom = "path") +
  #scale_x_continuous(breaks=c(20, 40, 80, 120), limits=c(20, 120))+
  theme_bw()+ labs(y=expression(R^2))+
  geom_hline(yintercept = 0)
#plots
ggsave("marginal_k_all_algos.pdf", plot=plots, device="pdf")

#Marginal effect of n
plotdat3$n <-factor(plotdat3$n, labels = c("40", "80", "160"))
plots <- ggplot(plotdat3, aes(x = n, y = R2, group=Algorithm, linetype = Algorithm, shape = Algorithm)) +
  stat_summary(fun.y = mean,
               geom = "point") +
  stat_summary(fun.y = mean,
               geom = "path") +
  #scale_x_continuous(breaks=c(40, 80, 160), limits=c(40, 160))+
  theme_bw()+ labs(y=expression(R^2))+
  geom_hline(yintercept = 0)
#plots
ggsave("marginal_n_all_algos.pdf", plot=plots, device="pdf")

#Marginal effect of tau2
plotdat3$tau2 <- factor(plotdat3$tau2, labels = c("0", "0.025", "0.05"))
plots <- ggplot(plotdat3, aes(x = tau2, y = R2, group=Algorithm, linetype = Algorithm, shape = Algorithm)) +
  stat_summary(fun.y = mean,
               geom = "point") +
  stat_summary(fun.y = mean,
               geom = "path") +
  #scale_x_continuous(breaks=c(0, .025, .05), limits=c(0, .05))+
  theme_bw()+ labs(y=expression(R^2), x=expression(tau^2))+
  geom_hline(yintercept = 0)
plots
ggsave("marginal_tau2_all_algos.pdf", plot=plots, device="pdf")



###################################################
#               Big R-squared plots               #
###################################################
includevars <- c("mf", "fx", "un", "ca")
plotdat <- data[ , .SD, .SDcols=c("k_studies", "mean_study_n", "es", "residual_heterogeneity", "moderators", "model", paste0(includevars, "_test_r2"))]
plotdat[, c(1:6):=lapply(.SD, factor), .SDcols=c(1:6)]

#Summarize several columns:
measure.vars <- names(plotdat)[-c(1:6)]
grouping.vars <- quote(list(k_studies,
                            mean_study_n,
                            es,
                            residual_heterogeneity,
                            moderators, model))

plotdat2<-plotdat[,lapply(.SD,mean, na.rm=TRUE),by=eval(grouping.vars), .SDcols=measure.vars]

id.vars=names(plotdat2)[1:6]
plotdat3<-melt.data.table(plotdat2, id.vars=id.vars)#measure.vars = measure.vars)
plotdat3[, "Algorithm" := substr(variable, 1,2)]
plotdat3[, variable := gsub("\\w+r2", "r2", variable)]

plotdat4<-dcast.data.table(plotdat3, k_studies + mean_study_n + es + residual_heterogeneity + moderators + model+Algorithm~variable, value.var="value")

names(plotdat4)[1:8]<-c("k", "n", "es", "tau2", "Moderators", "Model", "Algorithm", "R2")

categorical<-c("es", "tau2", "Moderators", "Model", "Algorithm")
plotdat4[, (categorical) := lapply(.SD, factor), .SDcols=categorical]
levels(plotdat4$Algorithm)
plotdat4$Algorithm <- ordered(plotdat4$Algorithm, levels = c("mf", "fx", "un", "ca"))
levels(plotdat4$Algorithm)<-c("MetaForest RE", "MetaForest FX", "MetaForest UN", "metaCART")

plotdat4$tau2 <- factor(plotdat4$tau2, labels = c("tau^2:~0", "tau^2:~0.04", "tau^2:~0.28"))
plotdat4$es <- factor(plotdat4$es, labels = c("beta:~0.2", "beta:~0.5", "beta:~0.8"))
plotdat4$n <- factor(plotdat4$n, labels = c("bar(n):~40", "bar(n):~80", "bar(n):~160"))
plotdat4$Model <- factor(plotdat4$Model, labels = c("a", "b", "c", "d", "e"))

plots<-lapply(levels(plotdat4$Model), function(x){
  tmp<-plotdat4[Model==x, -6]
  ggplot(tmp, aes(x = k, y = R2, linetype = Algorithm, shape = Algorithm, group=Algorithm)) +
    geom_point() +
    geom_path() +
    facet_grid(es + tau2 ~ Moderators+n, labeller = labeller(es = label_parsed, tau2 = label_parsed, n=label_parsed, Moderators=label_both))+
    theme_bw()+ labs(y=expression(R^2))+
    geom_hline(yintercept = 0)
})

for(i in 1:length(plots)){
  ggsave(paste0("Rsquare_plot_bw", i,".pdf"), plot=plots[[i]], device="pdf", width=297, height=210, units="mm")
}


##########################################
# Variable importance analyses and plots #
##########################################

#
#Merge files
#
res<-data.table(readRDS(list.files(pattern = "summarydata")))

n_chunks<-400
chunk_length<-405


res[,c(paste("mf_V", c(1:20), sep="")):=double() ]

for(i in 1:n_chunks){
  fileName<-paste0(c("metaforest_importance_", i, ".RData"), collapse="")
  if (file.exists(fileName)){
    imp.values<-readRDS(fileName)
    imp.values<-lapply(imp.values, function(x){c(x, rep(NA, 20-length(x)))})
    imp.values<-as.data.table(t(do.call(cbind, imp.values)))
    set(x=res, i=(1+((i-1)*chunk_length)):(i*chunk_length), j=c(paste("mf_V", c(1:20), sep="")), value=imp.values)
  }
}

#Metacart
res[,c(paste("ca_V", c(1:20), sep="")):=double() ]
for(i in 1:n_chunks){
  fileName<-paste0(c("cart_importance_", i, ".RData"), collapse="")
  if (file.exists(fileName)){
    imp.values<-readRDS(fileName)
    imp.values<-lapply(imp.values, function(x){
      full<-c(rep(0, 20))
      names(full)<-paste("V", c(1:20), sep="")
      full[names(x)]<-x
      full
    })
    imp.values<-as.data.table(t(do.call(cbind, imp.values)))
    set(x=data, i=(1+((i-1)*chunk_length)):(i*chunk_length), j=c(paste("ca_V", c(1:20), sep="")), value=imp.values)
  }
}

res[,nvars:=moderators]
for(nv in as.numeric(names(table(res$nvars)))[-3]){
  set(res, i=which(res$nvars==nv), j=c((min(grep("^ca_V", names(res)))+nv):max(grep("^ca_V", names(res)))), value=NA)
}


###
#Check which design factors are relevant for a correlated moderator
#
yvars<-list("mf_V1", "ca_V1")

anovas<-lapply(yvars, function(yvar){
  form<-paste(yvar, "(es + residual_heterogeneity + k_studies + mean_study_n + moderators + model) ^ 2", sep="~")
  thisaov<-aov(as.formula(form), data=res)
  thisetasq<-EtaSq(thisaov)[ , 2]
  list(thisaov, thisetasq)
})
etasqs<-data.frame(sapply(anovas, function(x){
  tempnums<-x[[2]]
  formatC(tempnums, 2, format="f")
}))

coef.table<-etasqs
names(coef.table)<-yvars

table.names<-row.names(coef.table)
table.names<-gsub("_studies", "", table.names)
table.names<-gsub("moderators", "M", table.names)
table.names<-gsub("residual_heterogeneity", "$\\tau^2$", table.names)
table.names<-gsub("es", "$\\beta$", table.names)
table.names<-gsub("mean_study_n", "$\\mean{n}$", table.names)
table.names<-gsub("model", "Model ", table.names)

row.names(coef.table)<-table.names
names(coef.table)<-c("MetaForest", "CA importance")

study1.coef.table<-coef.table
sink("output.txt", append = TRUE)
study1.coef.table
sink()

library(xtable)

print(xtable(study1.coef.table), file="table_study1_importance.tex",sanitize.text.function=function(x){x})


#Standardize importances
plotdat<-res

for(var in c("mf", "ca")){
  plotdat[, imp.sum := rowSums(abs(.SD), na.rm=TRUE), .SDcols=grep(paste0("^", var, "_V"), names(plotdat), value = TRUE)]

  for(thisvar in grep(paste0("^", var, "_V"), names(plotdat), value = TRUE)){
    plotdat[!(imp.sum == 0), (thisvar) := (.SD/imp.sum)*100, .SDcols = thisvar]
    #set(plotdat, i = !(imp.sum == 0), j=var, value=(plotdat[[var]]/plotdat[["imp.sum"]])*100)
  }
}

plotdat[, c("nvars", "imp.sum") := NULL]
#plotdat[k_studies==120 & mean_study_n == 160 & es == .8 & residual_heterogeneity == 0 & moderators == 5 & model == 4, ]
plotdat<-melt(plotdat, id.vars = names(plotdat)[c(1:7)])

plotdat[, Algorithm:=substr(variable, 1,2)]
library(stringr)
plotdat[, variable:= str_sub(variable, 5, -1)]

library(ggplot2)

names(plotdat)[2:9]<-c("k", "n", "es", "tau2", "Moderators", "Model", "Variable", "Importance")

plotdat[, c(2:8, 10) := lapply(.SD, factor), .SDcols=c(2:8, 10)]
#plotdat[k=="120" & n == "160" & es == "0.8" & tau2 == "0" & Moderators == "5" & Model == "4", ]
plotdat$es <- factor(plotdat$es, labels = c("beta:~0.2", "beta:~0.5", "beta:~0.8"))
plotdat$Model <- factor(plotdat$Model, labels = c("a", "b", "c", "d", "e"))
#plotdat[, Variable := factor(Variable, labels = c( "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20"))]
plotdat <- plotdat[Variable %in% c( "1", "2", "3", "4", "5"), ]
#plotdat <- plotdat[!(Algorithm == "cn"),]


plots<-lapply(levels(plotdat$Algorithm), function(x){
  tmp<-plotdat[Algorithm == x, -10]
  ggplot(tmp, aes(x = Variable, y = Importance)) +
    geom_boxplot(outlier.shape = NA) +
    facet_grid(Model + es ~ k, labeller =
                 labeller(Model = label_both, es = label_parsed, k=label_both)) +
    theme_bw()+geom_hline(yintercept = 0)
})

for(i in 1:length(plots)){
  ggsave(paste0("Importance_Plot_", levels(plotdat$Algorithm)[i],".pdf"), plot=plots[[i]], device="pdf", width=210, height=297, units="mm")
}


##########################################
# Tau2s                                  #
##########################################

#
#Merge files
#
res<-data.table(readRDS(list.files(pattern = "summarydata")))

n_chunks<-400
chunk_length<-405

###Read files
#Metaforest
res[,c("tau2"):=double() ]

for(i in 1:n_chunks){
  #i=279
  fileName<-paste0(c("res_tau2s_", i, ".RData"), collapse="")
  if (file.exists(fileName)){
    imp.values<-unlist(readRDS(fileName))
    set(x=res, i=(1+((i-1)*chunk_length)):(i*chunk_length), j=c("tau2"), value=imp.values)
  }
}


anovas <- aov(tau2 ~ (es + residual_heterogeneity + k_studies + mean_study_n + moderators + model) ^ 2, data=res)
etasqs <- EtaSq(anovas)[ , 2]
etasqs <- formatC(etasqs, 2, format="f")

coef.table<-etasqs
names(coef.table)<-yvars

plotdat <- res
names(plotdat)[2:8]<-c("k", "n", "es", "tau2", "Moderators", "Model", "tau2_est")

plotdat[, c(2:7) := lapply(.SD, factor), .SDcols=c(2:7)]

plotdat$es <- factor(plotdat$es, labels = c("beta:~0.2", "beta:~0.5", "beta:~0.8"))
plotdat$Model <- factor(plotdat$Model, labels = c("a", "b", "c", "d", "e"))
plotdat$tau2 <- factor(plotdat$tau2, labels = c("tau^2:~0", "tau^2:~0.04", "tau^2:~0.28"))
plotdat$n <- factor(plotdat$n, labels = c("bar(n):~40", "bar(n):~80", "bar(n):~160"))

plots<-lapply(levels(plotdat$Model), function(x){
  x <-levels(plotdat$Model)[1]
  tmp<-plotdat[Model == x, -7]
  ggplot(tmp, aes(x = tau2_est)) +
    geom_density()+
    geom_vline(aes(xintercept = as.numeric(gsub("tau\\^2\\:\\~", "", tau2))))+
    facet_grid(es + k ~ tau2, labeller =
                 labeller(es = label_parsed, tau2 = label_parsed, k=label_both)) +
    theme_bw()
})
for(i in 1:length(plots)){
  ggsave(paste0("Tau2_residual_Plot_", levels(plotdat$Model)[i],".pdf"), plot=plots[[i]], device="pdf", width=210, height=297, units="mm")
}f
