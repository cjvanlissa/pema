####################
# Data exploration #
####################
dependencies <- c('data.table', 'tidyverse', 'stringr', 'ggplot2')
lapply(dependencies, function(x) library(x, character.only = T))
dat <- as.data.table(readRDS("simulatie_2021/results/sim_results_2021-09-08.RData"))

# load in needed functions
source('./simulatie_2021/Analysis_functions.R') #created functions used for analysis

#some data preprocessing
newnames <- c('hs', 'lasso', 'mf_r', 'rma') #names of algorithms
conditions <- c("k_train", "mean_n", "es", "tau2","alpha_mod", "moderators", "model") #conditions that vary
lc <- length(conditions) #length of conditions, handy for further code
dat[, c(1:lc):=lapply(.SD, factor), .SDcols=c(1:lc)] #convert to factor
dat$model <- as.factor(dat$model)
levels(dat$model) <- c('exponential', 'linear', 'cubic', 'Twoway_interaction') #rename to more understandable levels


#there are duplicated columns which contain the same values, we can omit those
last <- function(x) { return( x[length(x)] ) } #function to obtain last element from vector
unique_cols <- which(!duplicated(colnames(dat)))
dat <- dat[ , unique_cols[1]:last(unique_cols)]

#Checking for missing values
#first we omit the columns that should contain missing values which are the 'sel' columns
lashs <- colnames(dat)[grep(pattern = "_(ci|hdi)_sel_\\d", colnames(dat))]
mfrma <- colnames(dat)[grep(pattern = "(mf|rma)_sel_\\d", colnames(dat))]

#create df without sel columns and subset it for rows that contain more than 1 missing
wo_sel_cols <- colnames(dat)[!colnames(dat) %in% c(lashs, mfrma)]
dat_no_sel <- dat[ , .SD, .SDcols=c(conditions, wo_sel_cols)] #subset data
MD <- dat_no_sel[rowSums(is.na(dat_no_sel)) > 0, ]
View(MD) # issue is with rma

#where is the issue?
apply(MD[, 1:lc], 2, table) #cubic and exponential model when M = 2,3,6

#Not all models are run with all levels for moderators
table(dat[dat$moderators == 7 | dat$moderators == 4,]$model) #note that moderators = 4 or 7 only exists for 2-way
table(dat[dat$moderators == 2| dat$moderators == 6,]$model)  #and when moderators = 2, 6 none exist for 2way
table(dat[dat$moderators == 3,]$model) #3 is ran in every model equally

#design factors we need for further code
grouping.vars <- quote(list(k_train,
                            mean_n,
                            es,
                            tau2,
                            alpha_mod,
                            moderators,
                            model))

##################################################
# Predictive performance analyses - MSE & R2     #
##################################################

#subsets the data based on particular performance criterion
#in this case test and train_r2
analyzedat_test <- dat[ , .SD, .SDcols=c(conditions, paste0(newnames, "_test_r2"))]
analyzedat_train <- dat[ , .SD, .SDcols=c(conditions, paste0(newnames, "_train_r2"))]

#delete the 20 cases for which RMA is missing
analyzedat_test <- na.omit(analyzedat_test)
analyzedat_train <- na.omit(analyzedat_train)

#convergence checks, traceplots
Traceplot(analyzedat_test$rma_test_r2, analyzedat_train$rma_train_r2, "R2", "RMA")
Traceplot(analyzedat_test$mf_r_test_r2, analyzedat_train$mf_r_train_r2, "R2", "Metaforest")
Traceplot(analyzedat_test$lasso_test_r2, analyzedat_train$lasso_train_r2, "R2", "Lasso")
Traceplot(analyzedat_test$hs_test_r2, analyzedat_train$hs_train_r2, "R2", "Horseshoe")

#convergence checks, densities
par(mfrow = c(2, 2))
plotdens(analyzedat_test)
plotdens(analyzedat_train)

#obtain descriptives for train and test R2
psych::describe(analyzedat_train[,(lc+1):ncol(analyzedat_train)])
psych::describe(analyzedat_test[,(lc+1):ncol(analyzedat_test)])


#the dependent variables (the test r2)
yvars <- paste0(newnames, "_test_r2")
#creates a list for every Anova with the results for the anova and the effect sizes for the conditions on the algorithms
anovas<-lapply(yvars, function(yvar){
  form<-paste(yvar, '~(', paste(unlist(conditions[-lc]), "+", collapse = ' '), conditions[lc], ") ^ 2") #the ^2 signifies that we want all possible effects up until interactions
  thisaov<-aov(as.formula(form), data=analyzedat_test) #change data according to dataframe you are using
  thisetasq<-EtaSq(thisaov)[ , 2]
  list(thisaov, thisetasq)
})

#creates dataframe with effect sizes for all conditions for all algorithms
etasqs<-data.frame(sapply(anovas, function(x){
  tempnums<-x[[2]]
  formatC(tempnums, 2, format="f")
}))

coef.table<-copy(etasqs)

# row.names(coef.table)<-table.names
names(coef.table)<-c("Horseshoe", "Lasso", "Metaforest" ,"RMA")
coef.table.sort <- coef.table[order(coef.table[,"Horseshoe"], decreasing = T), , drop = FALSE]


table.names<-row.names(coef.table.sort)
table.names<-gsub("moderators", "M", table.names) #this changes 'moderators' to M
table.names<-gsub("tau2", "$\\tau^2$", table.names) #nothing
table.names<-gsub("es", "$\\beta$", table.names) #changes 'es'to '$beta$ but I do not know why
table.names<-gsub("mean_n", "$\\mean{n}$", table.names) #nothing
table.names<-gsub("model", "Model ", table.names) #capitalizes 'model'
table.names<-gsub("k_train", "$\\kappa$ ", table.names)
table.names<-gsub("alpha_mod", "$\\alpha$ ", table.names)
row.names(coef.table.sort)<-table.names

write.csv(coef.table.sort, file =  "./simulatie_2021/Analysis Results/EtaSq_testR2_sort.csv")
etasq <- read.csv2(file =  "./simulatie_2021/Analysis Results/EtaSq_testR2.csv", sep = ",", row.names = 'X')

#subsets the data based on particular performance criterion
#in this case test and mse
analyzedat_test_mse <- dat[ , .SD, .SDcols=c(conditions, paste0(newnames, "_test_mse"))]
analyzedat_train_mse <- dat[ , .SD, .SDcols=c(conditions, paste0(newnames, "_train_mse"))]

analyzedat_test_mse <- na.omit(analyzedat_test_mse)
analyzedat_train_mse <- na.omit(analyzedat_train_mse)

#obtain descriptives for train and test mse
psych::describe(analyzedat_train_mse[,(lc+1):(lc+4)])
psych::describe(analyzedat_test_mse[,(lc+1):(lc+4)])



############
# R2 plots #
############
measure.vars <- names(analyzedat_test)[-c(1:lc)] #test_r2 values

#based on median
testR2_per_condition <- analyzedat_test[,lapply(.SD, median),by=eval(grouping.vars), .SDcols=measure.vars]

metrics_r2 <- names(testR2_per_condition)[grep("r2", colnames(testR2_per_condition))]

plot_interaction_median('es', 'model', testR2_per_condition, metrics_r2, "R2")
plot_marginal_median('model', testR2_per_condition, metrics_r2, "R2", pointsize = 5, linesize = 2)
plot_interaction_median('alpha_mod', 'model', testR2_per_condition, metrics_r2, "R2")
plot_marginal_median('tau2', testR2_per_condition, metrics_r2, "R2", 5, 2)
plot_marginal_median('mean_n', testR2_per_condition, metrics_r2, "R2", 5, 2)
plot_marginal_median('k_train', testR2_per_condition, metrics_r2, "R2", 5, 2)
plot_marginal_median('moderators', testR2_per_condition, metrics_r2, "R2", 5, 2)



#---End R2---------------------------------------------------------------------------

########
# tau2 #
########
#create tau2 difference variables which are the estimate tau2 - true tau2
dat <- dat %>%
  mutate(hs_diff_tau2 = hs_tau2 - as.numeric(as.character(tau2)),
         lasso_diff_tau2 = lasso_tau2 - as.numeric(as.character(tau2)),
         mf_r_diff_tau2 = mf_r_tau2 - as.numeric(as.character(tau2)),
         rma_diff_tau2 = rma_tau2 - as.numeric(as.character(tau2)))



#subset for conditions and the new tau2 diff variables
analyzedat_tau <- dat[ , .SD, .SDcols=c(conditions, paste0(newnames, "_diff_tau2"))]
#omit NA's which are still the 20 NA's from the rma algorithm
analyzedat_tau <- na.omit(analyzedat_tau)
#obtain descriptives
psych::describe(analyzedat_tau[,(lc+1):ncol(analyzedat_tau)])

#obtain names of diff variables and conditions
measure.vars <- names(analyzedat_tau)[-c(1:lc)]

#the dependent variables (the diff tau2)
yvars <- paste0(newnames, "_diff_tau2")
#creates a list for every Anova with the results for the anova and the effect size for the conditions on the algorithms
anovas<-lapply(yvars, function(yvar){
  form<-paste(yvar, '~(', paste(unlist(conditions[-lc]), "+", collapse = ' '), conditions[lc], ") ^ 2")
  thisaov<-aov(as.formula(form), data=analyzedat_tau)
  list(thisaov, thisetasq)
})

#creates dataframe with effect sizes for all conditions for all algorithms
etasqs<-data.frame(sapply(anovas, function(x){
  tempnums<-x[[2]]
  formatC(tempnums, 2, format="f")
}))

coef.table<-copy(etasqs)

# row.names(coef.table)<-table.names
names(coef.table)<-c("Horseshoe", "Lasso", "Metaforest" ,"RMA")
coef.table.sort <- coef.table[order(coef.table[,"Horseshoe"], decreasing = T), , drop = FALSE]


table.names<-row.names(coef.table.sort)
table.names<-gsub("moderators", "M", table.names) #this changes 'moderators' to M
table.names<-gsub("tau2", "τ2", table.names) #nothing
table.names<-gsub("es", "β", table.names) #changes 'es'to '$beta$ but I do not know why
table.names<-gsub("mean_n", "n", table.names) #nothing
table.names<-gsub("model", "Model ", table.names) #capitalizes 'model'
table.names<-gsub("k_train", "κ", table.names)
table.names<-gsub("alpha_mod", "α", table.names)
row.names(coef.table.sort)<-table.names

write.csv(coef.table.sort, file =  "./simulatie_2021/Analysis Results/EtaSq_tau2_sort.csv")

# plot based on median
tau2_per_condition <- analyzedat_tau[,lapply(.SD, median),by=eval(grouping.vars), .SDcols=measure.vars]
metrics_tau2 <- names(tau2_per_condition)[grep("diff_tau2", colnames(tau2_per_condition))]

plot_marginal_median('model', tau2_per_condition, metrics_tau2, "\U0394Tau2")
plot_interaction_median('es', 'model', tau2_per_condition, metrics_tau2, "\U0394Tau2")
plot_marginal_median('es',tau2_per_condition, metrics_tau2, "\U0394Tau2")
plot_interaction_median('alpha_mod', 'model', tau2_per_condition, metrics_tau2, "\U0394Tau2")
plot_marginal_median('tau2', tau2_per_condition, metrics_tau2, "\U0394Tau2")
plot_interaction_median('moderators', 'model', tau2_per_condition, metrics_tau2, "\U0394Tau2")
plot_interaction_median('k_train', 'model', tau2_per_condition, metrics_tau2, "\U0394Tau2")
plot_marginal_median('mean_n', tau2_per_condition, metrics_tau2, "\U0394Tau2")


#---End tau2-----------------------
rm(analyzedat_tau, analyzedat_test, analyzedat_train, analyzedat_test_mse, analyzedat_train_mse,
   dat_no_sel, MD, tau2_per_condition, testR2_per_condition) #to save up memory

######################
# Variable Selection #
######################

#lashs and mfrma have been defined in lines 27:31
varsel <- dat[ , .SD, .SDcols=c(conditions, lashs, mfrma)] #subset data for conditions and selection variables

#below we see that the variables selected by CI and HDI are exaclty the same for lasso and Horseshoe
hs_sel <- lashs[grep('hs', lashs)] #obtain vector with CI and HDI names for horseshoe
las_sel <- lashs[grep('las', lashs)] #same for lasso

samesel <- function(x,y){sum(x == y, na.rm = T) / length(x[!is.na(x)])} #give proportion of two vectors where x == y
for(i in 1:7){
  print(paste0('proportion CI and HDI chose same variables: ', samesel(dat[[hs_sel[i]]], dat[[hs_sel[i+7]]])))
  print(paste0('proportion CI and HDI chose same variables: ', samesel(dat[[las_sel[i]]], dat[[las_sel[i+7]]])))
} #they are all 1, meaning that HDI and CI columns are exactly the same


#this is done to cbind the dataframes created for the different algorithms later in the code
varsel <- varsel[order(varsel$moderators),]
varsel$identifier <- 1:nrow(varsel)

#apply sel_all_levels function to all patterns and store in list
#this returns a dataframe with true negatives and true positives for an algorithms
all_patterns <- c("hs_ci_sel_\\d","las_ci_sel_\\d", "rma_sel_\\d", "mf_sel_\\d")
MOAL <- lapply(all_patterns, sel_all_levels) #this takes a short while
names(MOAL) <- str_extract(all_patterns,"[:alpha:]+_[:alpha:]+") #give appropriate names

#create final df with all metrics for all algorithms
varsel_final <- MOAL[[1]] #set initial dataframe
for(i in 2:length(MOAL)){
  varsel_final <- suppressMessages(full_join(varsel_final, MOAL[[i]]))
}
rm(MOAL) #remove big object
varsel_final <- subset(varsel_final, select = -c(identifier)) #get rid of identifier variable
psych::describe(varsel_final[, (lc+1):ncol(varsel_final)]) #very high proportions TP and TN, except for metaforest TN


######################
# Plotting TN and TP #
######################
measure.vars <- colnames(varsel_final)[grep("(TP|TN)", colnames(varsel_final))]

#obtain mean TN and TP per condition
sel_per_condition <- varsel_final[,lapply(.SD, function(x){mean(x, na.rm = T)}),by=eval(grouping.vars), .SDcols=measure.vars]
sel_per_condition$model <- as.factor(sel_per_condition$model)

#plotting TN and TP
metrics_tn_names <- names(sel_per_condition)[grep("TN", colnames(sel_per_condition))]
metrics_tn <- c(metrics_tn_names[1], metrics_tn_names[2], metrics_tn_names[4], metrics_tn_names[3]) #had to swap here
#so that it works in plot functions

#first subset it for the model where true effect size != 0
sig_model <- sel_per_condition[sel_per_condition$es != 0,] #subset to model when es != 0
metrics_tp_names <- names(sig_model)[grep("TP", colnames(sig_model))]
metrics_tp <- c(metrics_tp_names[1], metrics_tp_names[2], metrics_tp_names[4], metrics_tp_names[3]) #here too


#plotsTN
plot_interaction_mean('es', 'model', sel_per_condition, metrics_tn, 'TN')
plot_interaction_mean('moderators', 'model', sel_per_condition, metrics_tn, 'TN')

#marginal plots TP
plot_marginal_mean('mean_n', sig_model, metrics_tp, 'TP')
plot_marginal_mean('tau2', sig_model, metrics_tp, 'TP')


#these seem to be noteworthy interactions TP
plot_interaction_mean('k_train', 'model', sig_model, metrics_tp, 'TP')
plot_interaction_mean('moderators', 'model', sig_model, metrics_tp, 'TP')



