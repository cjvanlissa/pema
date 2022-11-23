####################
# Data exploration #
####################
library(latex2exp)
library(cowplot)
library(data.table)
library(ggplot2)
library(tidySEM)
library(rcompanion)
# load in needed functions
source('./simulatie_2021/Analysis_functions.R') #created functions used for analysis

out <- list()

dat <- as.data.table(readRDS("simulatie_2021/results/sim_results_2021-09-08.RData"))
dat <- dat[dat$model %in% c("es * x[, 1]", "es * x[, 1] + es * (x[, 1] ^ 2) + es * (x[, 1] ^ 3)"),]
dat[, grep("^mf", names(dat)) := NULL]


#some data preprocessing
newnames <- c('hs', 'lasso', 'rma') #names of algorithms
conditions <- c("k_train", "mean_n", "es", "tau2","alpha_mod", "moderators", "model") #conditions that vary
lc <- length(conditions) #length of conditions, handy for further code
dat[, (conditions):=lapply(.SD, ordered), .SDcols=conditions] #convert to factor
levels(dat$model) <- c('linear', 'cubic') #rename to more understandable levels

#there are duplicated columns which contain the same values, we can omit those
unique_cols <- which(!duplicated(colnames(dat)))
dat <- dat[ , unique_cols[1]:tail(unique_cols, 1)]

#Checking for non converged models
dat_no_sel <- dat[ , .SD, .SDcols=c(conditions,
                                    grep("test_r2", names(dat), fixed = T, value = T))] #subset data
MD <- dat_no_sel[rowSums(is.na(dat_no_sel)) > 0, ]
out$missing <- nrow(MD)
#View(MD) # issue is with rma

#where is the issue?
#apply(MD[, 1:lc], 2, table) #cubic model

##################################################
# Predictive performance analyses - MSE & R2     #
##################################################

#subsets the data based on particular performance criterion
#in this case test and train_r2
analyzedat_test <- dat[ , .SD, .SDcols=c(conditions, paste0(newnames, "_test_r2"))]
analyzedat_test <- na.omit(analyzedat_test)

#the dependent variables (the test r2)
yvars <- paste0(newnames, "_test_r2")
#creates a list for every Anova with the results for the anova and the effect sizes for the conditions on the algorithms
anovas<-lapply(yvars, function(yvar){
  form<-paste(yvar, '~', paste(conditions, collapse = "+"))
  # form<-paste(yvar, '~(', paste(unlist(conditions[-lc]), "+", collapse = ' '), conditions[lc], ") ^ 2") #the ^2 signifies that we want all possible effects up until interactions
  thisaov<-aov(as.formula(form), data=analyzedat_test) #change data according to dataframe you are using
  thisetasq<-EtaSq(thisaov)[ , 2]
  list(thisaov, thisetasq)
})


# Anova for the difference ------------------------------------------------

comps <- expand.grid(yvars, yvars)
comps <- comps[!comps$Var1 == comps$Var2, ]
comps <- t(apply(comps, 1, sort))
comps <- comps[!duplicated(comps), ]
#creates a list for every Anova with the results for the anova and the effect sizes for the conditions on the algorithms
diffanovas <- sapply(1:nrow(comps), function(i){
  form<-as.formula(paste("r2", '~ algo * (', paste(conditions, collapse = "+"), ")")) #the ^2 signifies that we want all possible effects up until interactions
  # form<-as.formula(paste("r2", '~ algo * ((', paste(unlist(conditions[-lc]), "+", collapse = ' '), conditions[lc], ") ^ 2)")) #the ^2 signifies that we want all possible effects up until interactions
  tmp <- analyzedat_test
  tmp <- tmp[, .SD, .SDcols = c(names(tmp)[1:7], comps[i, , drop = TRUE])]
  names(tmp) <- gsub("^(.+?)_r2", "r2_\\1", names(tmp))
  tmp = melt(tmp, id.vars = names(tmp)[1:7],
             measure.vars = names(tmp)[8:9],
             variable.name = "algo",
             value.name = "r2")
  thisaov<-aov(form, data=tmp) #change data according to dataframe you are using
  thisetasq<-EtaSq(thisaov)[ , 2]
  thisetasq <- thisetasq[startsWith(names(thisetasq), "algo")]
  thisetasq
})
colnames(diffanovas) <- paste0(gsub("_.+$", "", toupper(comps[,1])), " vs. ", gsub("_.+$", "", toupper(comps[,2])))
diffanovas <- data.frame(condition = gsub("algo:", "", rownames(diffanovas), fixed = T), diffanovas)
diffanovas$condition <- trimws(diffanovas$condition)


out$difference <- diffanovas[1, ]
#creates dataframe with effect sizes for all conditions for all algorithms
etasqs<-data.frame(sapply(anovas, `[[`, 2))
etasqs$condition <- trimws(rownames(etasqs))
etasqs <- merge(etasqs, diffanovas, by = "condition", all.x = TRUE)

table_anova<-copy(etasqs)

# row.names(table_anova)<-table.names
names(table_anova)[2:4]<-c("HS", "LASSO", "RMA")
# table_anova <- table_anova[order(table_anova[,"HS"], decreasing = T), , drop = FALSE]


############
# R2 plots #
############

#based on median
r2bycond <- analyzedat_test[,lapply(.SD, median),by=conditions, .SDcols=yvars]
out$conditions <- nrow(r2bycond)

int_main <- sapply(names(r2bycond)[1:7], function(n){
  tmp <- r2bycond[, list(median(hs_test_r2), median(lasso_test_r2), median(rma_test_r2)), by=n]
  x <- sapply(as.data.frame(tmp)[2:4], interpret)
  if(length(unique(x)) == 1){
    return(x[1])
  } else{
    "other"
  }
  })
names(int_main) <- names(r2bycond)[1:7]

# ints_names <- grep(":", table_anova$condition, fixed = TRUE, value = TRUE)
# ints <- do.call(rbind, lapply(strsplit(ints_names, ":"), unlist))
# int_int <- mapply(function(n1, n2){
#   eff <- "spreading"
#   out <- r2bycond[, list(median(hs_test_r2), median(lasso_test_r2), median(rma_test_r2)), by=eval(parse(text = paste0("list(", n1, ",", n2, ")")))]
#   for(val in unique(out[[1]])){
#     tmp <- as.data.frame(out[out[[1]] == val, ])
#     tmp <- sapply(tmp[3:5], interpret)
#     if(!all(tmp == int_main[names(out)[2]])){
#       eff <- "crossover"
#     }
#   }
#   for(val in unique(out[[2]])){
#     tmp <- as.data.frame(out[out[[2]] == val, ])
#     tmp <- sapply(tmp[3:5], interpret)
#     if(!all(tmp == int_main[names(out)[1]])){
#       eff <- "crossover"
#     }
#   }
#   return(eff)
# }, n1 = ints[, 1], n2 = ints[,2])

# cond_name <- mapply(function(x, y){
#   out <- c(paste0(x, ":", y), paste0(y, ":", x))
#   out[which(out %in% diffanovas$condition)]
# }, x = ints$Var1, y = ints$Var2)

interpretations <- #rbind(
  data.frame(condition = names(int_main),
             Interpretation = int_main)#,
  # data.frame(condition = ints_names,
  #            Interpretation = int_int)
  # )
table_anova <- merge(table_anova, interpretations, by= "condition", all.x = TRUE)

table.names<-renamefactors(table_anova$condition)
table_anova$Factor<-table.names

saveRDS(table_anova, "paper/anova.RData")
write.csv(table_anova, file =  "paper/anova.csv", row.names = FALSE)

# plotthese <- table_anova$condition[table_anova$Interpretation == "crossover"]
# p <- lapply(1:length(plotthese), function(i){
#   x = plotthese[i]
#   namfct <- paste0(gsub("^.+?:", "", x), " (facets)")
#   namx <- gsub(":.+$", "", x)
#   namfct<-renamefactors(namfct)
#   namx<-renamefactors(namx)
#   namfct <- paste0(letters[i], ") ", namfct)
#   plot_interaction_median(gsub(":.+$", "", x), gsub("^.+?:", "", x), r2bycond, yvars, "R2") + labs(y = NULL, x = TeX(namx), title = TeX(namfct))
# })
#
# p[[1]] <- p[[1]] + theme(legend.position = c(.76, .22), legend.title = element_blank())
# #p[[1]] <- p[[1]] + theme(legend.position = c(.16, .22), legend.title = element_blank())
# plot_other <- do.call(plot_grid, c(p, list(labels = NULL, ncol = 2))) # "auto"
# #ggsave("paper/other_effects.png", plot = plot_other, device = "png", height = 210, width =  297, units = "mm", scale = 2.5)
# ggsave("paper/other_effects.png", plot = plot_other, device = "png", width = 210, height =  297, units = "mm", scale = 2.5)



# Plot greatest difference ------------------------------------------------
plotthese <- conditions[order(table_anova$HS.vs..RMA[match(conditions, table_anova$condition)], decreasing = T)]
p <- lapply(1:length(plotthese), function(i){
  x = plotthese[i]
  namx<-renamefactors(x)
  plot_marginal_median(x, r2bycond, yvars, "R2") + labs(y = NULL, x = TeX(namx), title = paste0(letters[i], ")")) + theme(legend.position = "none")
})

p[[1]] <- p[[1]] + theme(legend.position = c(.7, .22), legend.title = element_blank())
plot_other <- do.call(plot_grid, c(p, list(labels = NULL, ncol = 2))) # "auto"
#ggsave("paper/other_effects.png", plot = plot_other, device = "png", height = 210, width =  297, units = "mm", scale = 2.5)
ggsave("paper/r2.png", plot = plot_other, device = "png", width = 210, height =  297, units = "mm", scale = 2.5)

# How often each highest?
df_hi <- dat
which_highest <- table(apply(df_hi[!es == "0", .SD, .SDcols = yvars], 1, which.max))
names(which_highest) <- yvars
out$which_highest <- which_highest

# which_highest <- yvars[apply(df_hi[, .SD, .SDcols = yvars], 1, which.max)]
# ishs <- which_highest == "hs_test_r2"
# df_hi[, "ishs" := ishs]
# df_hi <- df_hi[, .SD, .SDcols = c(conditions, "ishs")]
# df_hi <- df_hi[,lapply(.SD, sum),by=conditions, .SDcols="ishs"]

r2mean <- rbind(
  overall = colMeans(r2bycond[, .SD, .SDcols = yvars]),
  notnull = colMeans(r2bycond[!es == "0", .SD, .SDcols = yvars]),
  null = colMeans(r2bycond[es == "0", .SD, .SDcols = yvars]))
colnames(r2mean) <- paste("Mean R2", c("HS", "LASSO", "RMA"))
r2sd <- rbind(
  overall = apply(r2bycond[, .SD, .SDcols = yvars], 2, sd),
  notnull = apply(r2bycond[!es == "0", .SD, .SDcols = yvars], 2, sd),
  null = apply(r2bycond[es == "0", .SD, .SDcols = yvars], 2, sd))
colnames(r2sd) <- paste("SD R2", c("HS", "LASSO", "RMA"))
tmp <- apply(r2bycond[, .SD, .SDcols = yvars], 2, quantile, probs = c(.025, .975))
tmp2 <- apply(r2bycond[!es == "0", .SD, .SDcols = yvars], 2, quantile, probs = c(.025, .975))
tmp3 <- apply(r2bycond[es == "0", .SD, .SDcols = yvars], 2, quantile, probs = c(.025, .975))
cis <- rbind(
  overall = conf_int(lb = tmp[1, ], ub = tmp[2, ]),
  notnull = conf_int(lb = tmp2[1, ], ub = tmp2[2, ]),
  null = conf_int(lb = tmp3[1, ], ub = tmp3[2, ]))
colnames(cis) <- paste("CI ", c("HS", "LASSO", "RMA"))
tab <- cbind(r2mean, cis, r2sd)
tab <- data.frame(ES = c("Overall", "ES = 0", "ES != 0"), tab, check.names = F)
write.csv(tab, "paper/r2.csv", row.names = F)


#---End R2---------------------------------------------------------------------------




######################
# Variable Selection #
######################

#below we see that the variables selected by CI and HDI are exaclty the same for lasso and Horseshoe
table(dat$hs_ci_sel_1, dat$hs_hdi_sel_1)
table(dat$hs_ci_sel_2, dat$hs_hdi_sel_2)

varsel <- dat[!es == "0", .SD, .SDcols = c(conditions, grep("(rma|hs_ci|las_ci)_sel_[12]", names(dat), value = TRUE))]
varsel <- na.omit(varsel)

selvars <- grep("(rma|hs_ci|las_ci)_sel_1", names(dat), value = TRUE)
notselvars <- grep("(rma|hs_ci|las_ci)_sel_2", names(dat), value = TRUE)


tntp <- varsel[,lapply(.SD, sum),by=eval(conditions), .SDcols=grep("(rma|hs_ci|las_ci)_sel_[12]", names(dat), value = TRUE)]
table(varsel$mean_n)
sel_prop <- unlist(tntp[, lapply(.SD, sum, na.rm = T), .SDcols=grep("(rma|hs_ci|las_ci)_sel_[12]", names(dat), value = TRUE)] / nrow(varsel))
sel_prop <- matrix(sel_prop, nrow = 2)
sel_prop[2, ] <- 1-sel_prop[2, ]

out$selection <- sel_prop


tab_select <- sapply(conditions, function(thiscond){
  tab <- as.data.frame(varsel[, lapply(.SD,sum), by=thiscond, .SDcols = selvars][, -1])
  tots <- table(varsel[[thiscond]])
  tots <- tots[!tots == 0]
  c(sapply(tab, function(x){ cramerV(rbind(x, tots-x)) }),
  mapply(function(i, j){
    cramerV(as.matrix(tab[c(i,j)]))
  }, i = 1:3, j = c(2,3,2)))
})
tab_select <- data.frame(Factor = renamefactors(colnames(tab_select)), t(tab_select))
names(tab_select)[-1] <- paste0("$P_{", names(table_anova)[2:7], "}$")

int_main <- sapply(rownames(tab_select), function(n){
  tmp <- varsel[,lapply(.SD, sum),by=n, .SDcols=selvars]
  x <- sapply(as.data.frame(tmp)[2:4], interpret)
  if(length(unique(x)) == 1){
    return(x[1])
  } else{
    "other"
  }
})
names(int_main) <- rownames(tab_select)
tab_select$Interpretation <- int_main
tab_select2 <- sapply(conditions, function(thiscond){
  tab <- as.data.frame(varsel[, lapply(.SD,sum), by=thiscond, .SDcols = notselvars][, -1])
  tots <- table(varsel[[thiscond]])
  tots <- tots[!tots == 0]
  c(sapply(tab, function(x){ cramerV(rbind(x, tots-x)) }),
    mapply(function(i, j){
      cramerV(as.matrix(tab[c(i,j)]))
    }, i = 1:3, j = c(2,3,2)))
})
tab_select2 <- data.frame(Factor = renamefactors(colnames(tab_select2)), t(tab_select2))
names(tab_select2)[-1] <- paste0("$N_{", names(table_anova)[2:7], "}$")

int_main <- sapply(rownames(tab_select2), function(n){
  tmp <- varsel[,lapply(.SD, sum),by=n, .SDcols=notselvars]
  x <- sapply(as.data.frame(tmp)[2:4], interpret)
  if(length(unique(x)) == 1){
    return(x[1])
  } else{
    "other"
  }
})
names(int_main) <- rownames(tab_select2)
tab_select2$Interpretation <- int_main

saveRDS(tab_select, "paper/selected.RData")
write.csv(tab_select, file =  "paper/selected.csv", row.names = FALSE)
saveRDS(tab_select2, "paper/notselected.RData")
write.csv(tab_select2, file =  "paper/notselected.csv", row.names = FALSE)


p <- lapply(1:length(conditions), function(i){
  x = plotthese[i]
  namx<-renamefactors(x)
  tots <- table(varsel[[x]])
  tots <- tots[!tots==0]
  df_plot <- as.data.frame(varsel[,lapply(.SD, sum),by=x, .SDcols=selvars])
  names(df_plot)[-1] <- paste0("Sensitivity.", c("HS", "LASSO", "RMA"))
  df_plot[-1] <- lapply(df_plot[-1], `/`, tots)
  df_plot <- reshape(df_plot, direction = "long", varying = names(df_plot)[-1], timevar = "Algorithm")
  ggplot(data = df_plot, aes_string(x = x, y = "Sensitivity",
                                    linetype = "Algorithm",
                                    group = "Algorithm",
                                    shape = "Algorithm")) +
    geom_line(size = 1.5) +
    geom_point(size = 5) +
    theme_bw(base_size = 25) +
    labs(y = NULL, x = TeX(namx), title = paste0(letters[i], ")")) + theme(legend.position = "none")+
    scale_y_continuous(limits = c(.72, 1))
})

p[[1]] <- p[[1]] + theme(legend.position = c(.7, .22), legend.title = element_blank())
plot_other <- do.call(plot_grid, c(p, list(labels = NULL, ncol = 2))) # "auto"
#ggsave("paper/other_effects.png", plot = plot_other, device = "png", height = 210, width =  297, units = "mm", scale = 2.5)
ggsave("paper/sensitivity.png", plot = plot_other, device = "png", width = 210, height =  297, units = "mm", scale = 2.5)

p <- lapply(1:length(conditions), function(i){
  x = plotthese[i]
  namx<-renamefactors(x)
  tots <- table(varsel[[x]])
  tots <- tots[!tots==0]
  df_plot <- as.data.frame(varsel[,lapply(.SD, sum),by=x, .SDcols=notselvars])
  df_plot[-1] <- lapply(df_plot[-1], function(x){ 1-(x/tots)})
  names(df_plot)[-1] <- paste0("Sensitivity.", c("HS", "LASSO", "RMA"))
  df_plot <- reshape(df_plot, direction = "long", varying = names(df_plot)[-1], timevar = "Algorithm")
  ggplot(data = df_plot, aes_string(x = x, y = "Sensitivity",
                                    linetype = "Algorithm",
                                    group = "Algorithm",
                                    shape = "Algorithm")) +
    geom_line(size = 1.5) +
    geom_point(size = 5) +
    theme_bw(base_size = 25) +
    labs(y = NULL, x = TeX(namx), title = paste0(letters[i], ")")) + theme(legend.position = "none") +
    scale_y_continuous(limits = c(.9, 1))
})

p[[1]] <- p[[1]] + theme(legend.position = c(.7, .22), legend.title = element_blank())
plot_other <- do.call(plot_grid, c(p, list(labels = NULL, ncol = 2))) # "auto"
#ggsave("paper/other_effects.png", plot = plot_other, device = "png", height = 210, width =  297, units = "mm", scale = 2.5)
ggsave("paper/specificity.png", plot = plot_other, device = "png", width = 210, height =  297, units = "mm", scale = 2.5)



########
# tau2 #
########
# df_cond <- dat[, .SD, .SDcols = conditions]
# df_cond[, names(df_cond)[1:6] := lapply(.SD, function(i)as.numeric(as.character(i))), .SDcols = names(df_cond)[1:6]]
# df_cond[, "moderators" := 1]
# df_cond[, "tau2" := 0]
# df_cond[, "k_train" := 1000]
# df_cond <- df_cond[!duplicated(df_cond), ]
# levels(df_cond$model) <- c("es * x[, 1]", "es * x[, 1] + es * (x[, 1] ^ 2) + es * (x[, 1] ^ 3)")
# source("./simulatie_2021/functions for simulation/simulate_smd.R")
# library(sn)
# library(metafor)
# empiricalvalues2 <- sapply(1:nrow(df_cond), function(i){
#   thisrow <- df_cond[i, ]
#   data <- do.call(simulate_smd, as.list(thisrow))
#   data <- data$training
#   res <- rma(yi = data$yi,
#              vi = data$vi,
#              mods = ~x, data = data)
#   c(res$b[2,1], res$tau2)
# #   data$study <- 1:10000
# #   tryCatch({tmp <- lme(yi ~ x, random = ~ 1 | study, weights = varFixed(~ vi), control=lmeControl(sigma = 1), data=data)
# #   c(tmp$coefficients$fixed[2],
# #   attr(summary(tmp$modelStruct$reStruct$study), "stdDev")^2)}, error = function(e){c(NA, NA)})
# })
#df_cond <- cbind(df_cond, t(empiricalvalues2))
measure.vars <- grep("_tau2$", names(dat), value = TRUE)
tau2 <- dat[, .SD, .SDcols = c(conditions, measure.vars)]
tau2[, (measure.vars[1:2]) := .SD^2, .SDcols = measure.vars[1:2]]
tau2[, (grep("_tau2$", names(tau2), value = TRUE)) := .SD - as.numeric(as.character(tau2)), .SDcols = grep("_tau2$", names(tau2), value = TRUE)]
tau2 <- na.omit(tau2)
out$tau2_bias <- tau2[,lapply(.SD, mean), .SDcols=measure.vars]
out$tau2_variance <- tau2[,lapply(.SD, var), .SDcols=measure.vars]
tau2_per_condition <- tau2[,lapply(.SD, median),by=eval(conditions), .SDcols=measure.vars]


yvars <- paste0(c("hs", "lasso", "rma"), "_tau2")
anovas<-lapply(yvars, function(yvar){
  form<-paste(yvar, '~', paste(conditions, collapse = " + "))
  #form<-paste(yvar, '~(', paste(unlist(conditions[-lc]), "+", collapse = ' '), conditions[lc], ") ^ 2") #the ^2 signifies that we want all possible effects up until interactions
  thisaov<-aov(as.formula(form), data=tau2)
  thisetasq<-EtaSq(thisaov)[ , 2]
  list(thisaov, thisetasq)
})


# Anova for the difference ------------------------------------------------

comps <- expand.grid(yvars, yvars)
comps <- comps[!comps$Var1 == comps$Var2, ]
comps <- t(apply(comps, 1, sort))
comps <- comps[!duplicated(comps), ]
diffanovas <- sapply(1:nrow(comps), function(i){
  form<-as.formula(paste("tau2hat", '~ algo * (', paste(conditions, collapse = "+"),")")) #the ^2 signifies that we want all possible effects up until interactions
  tmp <- tau2
  tmp <- tmp[, .SD, .SDcols = c(names(tmp)[1:7], comps[i, , drop = TRUE])]
  names(tmp) <- gsub("^(.+?)_tau2", "tau2hat_\\1", names(tmp))
  tmp = melt(tmp, id.vars = names(tmp)[1:7],
             measure.vars = names(tmp)[8:9],
             variable.name = "algo",
             value.name = "tau2hat")
  thisaov<-aov(form, data=tmp) #change data according to dataframe you are using
  thisetasq<-EtaSq(thisaov)[ , 2]
  thisetasq <- thisetasq[startsWith(names(thisetasq), "algo")]
  thisetasq
})
colnames(diffanovas) <- paste0(gsub("_.+$", "", toupper(comps[,1])), " vs. ", gsub("_.+$", "", toupper(comps[,2])))
diffanovas <- data.frame(condition = gsub("algo:", "", rownames(diffanovas), fixed = T), diffanovas)
diffanovas$condition <- trimws(diffanovas$condition)

out$tau2difference <- diffanovas[1, ]
#creates dataframe with effect sizes for all conditions for all algorithms
etasqs<-data.frame(sapply(anovas, `[[`, 2))
etasqs$condition <- trimws(rownames(etasqs))
etasqs <- merge(etasqs, diffanovas, by = "condition", all.x = TRUE)

table_tau2<-copy(etasqs)
table_tau2$Factor<- renamefactors(table_tau2$condition)
names(table_tau2)[2:4]<-c("HS", "LASSO", "RMA")

saveRDS(table_tau2, "paper/table_tau2.RData")
write.csv(table_tau2, file =  "paper/table_tau2.csv", row.names = FALSE)

# Variance of tau2 --------------------------------------------------------

tau_var <- tau2
tau_var <- tau_var[ , lapply(.SD, var), by=conditions, .SDcols = measure.vars]

yvars <- measure.vars
anovas<-lapply(yvars, function(yvar){
  form<-paste(yvar, '~', paste(conditions, collapse = " + "))
  #form<-paste(yvar, '~(', paste(unlist(conditions[-lc]), "+", collapse = ' '), conditions[lc], ") ^ 2") #the ^2 signifies that we want all possible effects up until interactions
  thisaov<-aov(as.formula(form), data=tau_var)
  thisetasq<-EtaSq(thisaov)[ , 2]
  list(thisaov, thisetasq)
})


# Anova for the difference ------------------------------------------------

comps <- expand.grid(yvars, yvars)
comps <- comps[!comps$Var1 == comps$Var2, ]
comps <- t(apply(comps, 1, sort))
comps <- comps[!duplicated(comps), ]
diffanovas <- sapply(1:nrow(comps), function(i){
  form<-as.formula(paste("tauhat", '~ algo * (', paste(conditions, collapse = "+"),")")) #the ^2 signifies that we want all possible effects up until interactions
  tmp <- tau_var
  tmp <- tmp[, .SD, .SDcols = c(conditions, comps[i, , drop = TRUE])]
  names(tmp) <- gsub("^(.+?)_tau1", "tauhat_\\1", names(tmp))
  tmp = melt(tmp, id.vars = names(tmp)[1:7],
             measure.vars = names(tmp)[8:9],
             variable.name = "algo",
             value.name = "tauhat")
  thisaov<-aov(form, data=tmp) #change data according to dataframe you are using
  thisetasq<-EtaSq(thisaov)[ , 2]
  thisetasq <- thisetasq[startsWith(names(thisetasq), "algo")]
  thisetasq
})
colnames(diffanovas) <- paste0(gsub("_.+$", "", toupper(comps[,1])), " vs. ", gsub("_.+$", "", toupper(comps[,2])))
diffanovas <- data.frame(condition = gsub("algo:", "", rownames(diffanovas), fixed = T), diffanovas)
diffanovas$condition <- trimws(diffanovas$condition)

#creates dataframe with effect sizes for all conditions for all algorithms
etasqs<-data.frame(sapply(anovas, `[[`, 2))
etasqs$condition <- trimws(rownames(etasqs))
etasqs <- merge(etasqs, diffanovas, by = "condition", all.x = TRUE)

table_tau<-copy(etasqs)
table_tau$Factor <- renamefactors(table_tau$condition)
names(table_tau)[2:4]<-c("HS", "LASSO", "RMA")

saveRDS(table_tau, "paper/table_tau_var.RData")
write.csv(table_tau, file =  "paper/Supplemental_table_S2_tau_variance.csv", row.names = FALSE)
usethis::use_git_ignore("!paper/Supplemental_table_S2_tau_variance.csv")


# Betas -------------------------------------------------------------------

measure.vars <- grep("_beta1$", names(dat), value = TRUE)
beta <- dat[, .SD, .SDcols = c(conditions, measure.vars)]
beta <- na.omit(beta)
beta[, centerbeta := c(0, .2, .5, .8)[as.integer(es)]]
beta[model == "cubic", centerbeta := c(0, 1.1, 2.6, 4)[as.integer(es)]]
beta[, (measure.vars) := .SD - centerbeta, .SDcols = measure.vars]

out$beta_bias <- beta[,lapply(.SD, mean), .SDcols=measure.vars]
out$beta_variance <- beta[,lapply(.SD, var), .SDcols=measure.vars]

yvars <- paste0(c("hs", "lasso", "rma"), "_beta1")
anovas<-lapply(yvars, function(yvar){
  form<-paste(yvar, '~', paste(conditions, collapse = " + "))
  #form<-paste(yvar, '~(', paste(unlist(conditions[-lc]), "+", collapse = ' '), conditions[lc], ") ^ 2") #the ^2 signifies that we want all possible effects up until interactions
  thisaov<-aov(as.formula(form), data=beta)
  thisetasq<-EtaSq(thisaov)[ , 2]
  list(thisaov, thisetasq)
})


# Anova for the difference ------------------------------------------------

comps <- expand.grid(yvars, yvars)
comps <- comps[!comps$Var1 == comps$Var2, ]
comps <- t(apply(comps, 1, sort))
comps <- comps[!duplicated(comps), ]
diffanovas <- sapply(1:nrow(comps), function(i){
  form<-as.formula(paste("betahat", '~ algo * (', paste(conditions, collapse = "+"),")")) #the ^2 signifies that we want all possible effects up until interactions
  tmp <- beta
  tmp <- tmp[, .SD, .SDcols = c(conditions, comps[i, , drop = TRUE])]
  names(tmp) <- gsub("^(.+?)_beta1", "betahat_\\1", names(tmp))
  tmp = melt(tmp, id.vars = names(tmp)[1:7],
             measure.vars = names(tmp)[8:9],
             variable.name = "algo",
             value.name = "betahat")
  thisaov<-aov(form, data=tmp) #change data according to dataframe you are using
  thisetasq<-EtaSq(thisaov)[ , 2]
  thisetasq <- thisetasq[startsWith(names(thisetasq), "algo")]
  thisetasq
})
colnames(diffanovas) <- paste0(gsub("_.+$", "", toupper(comps[,1])), " vs. ", gsub("_.+$", "", toupper(comps[,2])))
diffanovas <- data.frame(condition = gsub("algo:", "", rownames(diffanovas), fixed = T), diffanovas)
diffanovas$condition <- trimws(diffanovas$condition)

out$betadiff <- diffanovas[1, ]
#creates dataframe with effect sizes for all conditions for all algorithms
etasqs<-data.frame(sapply(anovas, `[[`, 2))
etasqs$condition <- trimws(rownames(etasqs))
etasqs <- merge(etasqs, diffanovas, by = "condition", all.x = TRUE)

table_beta<-copy(etasqs)
table_beta$Factor <- renamefactors(table_beta$condition)
names(table_beta)[2:4]<-c("HS", "LASSO", "RMA")

saveRDS(table_beta, "paper/table_beta.RData")
write.csv(table_beta, file =  "paper/table_beta.csv", row.names = FALSE)


# Variance of beta --------------------------------------------------------

beta_var <- beta
beta_var <- beta_var[ , lapply(.SD, var), by=conditions, .SDcols = measure.vars]

yvars <- measure.vars
anovas<-lapply(yvars, function(yvar){
  form<-paste(yvar, '~', paste(conditions, collapse = " + "))
  #form<-paste(yvar, '~(', paste(unlist(conditions[-lc]), "+", collapse = ' '), conditions[lc], ") ^ 2") #the ^2 signifies that we want all possible effects up until interactions
  thisaov<-aov(as.formula(form), data=beta_var)
  thisetasq<-EtaSq(thisaov)[ , 2]
  list(thisaov, thisetasq)
})


# Anova for the difference ------------------------------------------------

comps <- expand.grid(yvars, yvars)
comps <- comps[!comps$Var1 == comps$Var2, ]
comps <- t(apply(comps, 1, sort))
comps <- comps[!duplicated(comps), ]
diffanovas <- sapply(1:nrow(comps), function(i){
  form<-as.formula(paste("betahat", '~ algo * (', paste(conditions, collapse = "+"),")")) #the ^2 signifies that we want all possible effects up until interactions
  tmp <- beta_var
  tmp <- tmp[, .SD, .SDcols = c(conditions, comps[i, , drop = TRUE])]
  names(tmp) <- gsub("^(.+?)_beta1", "betahat_\\1", names(tmp))
  tmp = melt(tmp, id.vars = names(tmp)[1:7],
             measure.vars = names(tmp)[8:9],
             variable.name = "algo",
             value.name = "betahat")
  thisaov<-aov(form, data=tmp) #change data according to dataframe you are using
  thisetasq<-EtaSq(thisaov)[ , 2]
  thisetasq <- thisetasq[startsWith(names(thisetasq), "algo")]
  thisetasq
})
colnames(diffanovas) <- paste0(gsub("_.+$", "", toupper(comps[,1])), " vs. ", gsub("_.+$", "", toupper(comps[,2])))
diffanovas <- data.frame(condition = gsub("algo:", "", rownames(diffanovas), fixed = T), diffanovas)
diffanovas$condition <- trimws(diffanovas$condition)

#creates dataframe with effect sizes for all conditions for all algorithms
etasqs<-data.frame(sapply(anovas, `[[`, 2))
etasqs$condition <- trimws(rownames(etasqs))
etasqs <- merge(etasqs, diffanovas, by = "condition", all.x = TRUE)

table_beta<-copy(etasqs)
table_beta$Factor <- renamefactors(table_beta$condition)
names(table_beta)[2:4]<-c("HS", "LASSO", "RMA")

saveRDS(table_beta, "paper/table_beta_var.RData")
write.csv(table_beta, file =  "paper/Supplemental_table_S1_beta_variance.csv", row.names = FALSE)
usethis::use_git_ignore("!paper/Supplemental_table_S1_beta_variance.csv")


# Save all output ---------------------------------------------------------



saveRDS(out, "paper/output.RData")
