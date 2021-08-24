
##################################################
# Predictive performance analyses - Descriptives #
##################################################
library(data.table)
dat <- as.data.table(readRDS(r"(C:\Users\e_lib\OneDrive\Documents\Caspar Repositories\simulatie_2021\sim_results_2021-08-11.RData)"))
newnames <- c('hs', 'lasso', 'mf_r', 'rma') #names of models
conditions <- c("k_train", "mean_n", "es", "tau2","alpha_mod", "moderators", "model") #conditions that vary
lc <- length(conditions) #length of conditions, handy for further code

#subsets the data based on particular performance criterion
#in this case train_r2 and test_r2
analyzedat_train <- dat[ , .SD, .SDcols=c(conditions, paste0(newnames, "_train_r2"))]
analyzedat_test <- dat[ , .SD, .SDcols=c(conditions, paste0(newnames, "_test_r2"))]


#converts the condition variables to factors
analyzedat_train[, c(1:lc):=lapply(.SD, factor), .SDcols=c(1:7)]
analyzedat_test[, c(1:lc):=lapply(.SD, factor), .SDcols=c(1:7)]


setwd("./simulatie_2021")
#sink("output.txt") #makes connection with a .txt file, use it later
#unlink('output.txt') #unlinks the connection and deletes .txt file.


#there are missing values in the data for train
library(mice) #load in library

#show where missings are and prodivde summary of the dataset to check for patterns
misdat_pat <- function(df){
  varmis <- (apply(df, 2, function(x){sum(is.na(x))})) #shows where missings are
  cases_mis <- which(is.na(df)) #which cases are missing
  missdat <- df[cases_rma_mis,] #subset data for missing cases
  varmis_pat <- summary(missdat)
  mispat <- md.pattern(missdat, plot = T) #shows where missings are and how many
  return(list(variables_missing = varmis, pattern = varmis_pat))
}

misdat_pat(analyzedat_train) #17 missing in rma_train r2
misdat_pat(analyzedat_test) #the same for the test data


#Since only 17 cases out of 4000000 are missing, I do not think it will have a big impact on further analyses
#Although it might be interesting to see why the cases are missing, for now I will continue with the observed data
analyzedat_train <- na.omit(analyzedat_train)
analyzedat_test <- na.omit(analyzedat_test)

#obtain mean and standard deviations
mean_sd <- function(df){
  tmp<-rbind(round(df[, lapply(.SD, mean), .SDcols=c((lc+1):ncol(df))], 2),
             round(df[, lapply(.SD, sd), .SDcols=c((lc+1):ncol(df))], 2))
  tmp<-data.frame(names(tmp), t(tmp))
  apply(tmp, 1, function(x){paste(x, collapse=", ")})
  colnames(tmp) <- c('alg', 'mean', 'sd')
  return(tmp)
}
tmp_train <- mean_sd(analyzedat_train)
View(tmp_train)

tmp_test <- mean_sd(analyzedat_test)
View(tmp_test)

#this functions plots densities for the chains and provides descriptives
pre_analysis <- function(subdat = analyzedat, lb = -Inf, ub = Inf){ #lb = lowerbound, ub = upperbound for the boundaries of impossible values
  #plot densities for chains in one window
  par(mfrow = c(2, 2)) #we have 4 algorithms, so 4 plots will be produced
  subdat <- as.data.frame(subdat) #convert to dataframe
  for(column in (lc+1):ncol(subdat)){ #plot densities for all algorithms
    col_obs <- na.omit(subdat[,column]) #there are sometimes missings so omit them
    prop.impvals <- sum(col_obs < lb | col_obs > ub) / length(col_obs) #proportion of theoretically impossible values
    plot(density(col_obs),
         main = paste0('histogram of ', colnames(subdat)[column]),
         xlab = paste0('proportion impossible values: ', round(prop.impvals, 3)),
         xlim = c(-30,1)
         )
  }
  print(psych::describe(subdat[,(lc+1):ncol(subdat)])) #print descriptives of chains
}

pre_analysis(analyzedat_train, lb = 0 , ub = 1) #chains fo rma and mf look fine, not so much for hs and lasso.
#also, both hs and lasso have ~60% of r2 value that are impossible to obtain.
pre_analysis(analyzedat_test, lb = 0 , ub = 1) #both hs and lasso have ~55% of r2 value that should be impossible to obtain.
#although both mf and rma also contain impossible values.

#note that no algorithms produces r2 above 1, but the lowest value for rma is -0.33, so lets take that as a reference
sum(analyzedat_test$hs_test_r2 < -0.02) / nrow(analyzedat_test) #31.2%
#more severe outliers
sum(analyzedat_test$hs_test_r2 < -1) / nrow(analyzedat_test) #still 23%

#Which conditions show most of these problems?
#create 1 if r2 < -0.05
analyzedat_test$hs_fault <- ifelse(analyzedat_test$hs_test_r2 < -0.05, 1, 0)
analyzedat_test$lasso_fault <- ifelse(analyzedat_test$lasso_test_r2 < -0.05, 1, 0)


# run logistic regression on lasso, takes a while
# logreglasso <- glm(paste("lasso_fault", '~(', paste(unlist(conditions[-lc]), "+", collapse = ' '), conditions[lc], ") ^ 2")
# , data = analyzedat_test, family = 'binomial')

library(broom)
dfloglasso <- tidy(logreglasso)
dfloglasso <- na.omit(dfloglasso)
dfloglasso$positive <- ifelse(dfloglasso$estimate > 0, 1, 0)
dfloglasso$sig <- ifelse(dfloglasso$p.value < 0.001, 1, 0)
relevant_cond <- dfloglasso[dfloglasso$positive == 1 & dfloglasso$sig == 1,]
saveRDS(dfloglasso, file = 'logistic_reg_lasso.RData')
### ---- subset for positive coefficients to see whih wan increase the odds of R > 0.05
### --- also subset for when pvalue < 0.001



#subsetting data for when r2 < 0.05
conditions_lasso <- analyzedat_test[analyzedat_test$lasso_fault == 1, ]
conditions_lasso <- as.data.frame(conditions_lasso)


########################################
#Predictive performance on testing data#
########################################

analyzedat_test <- dat[ , .SD, .SDcols=c(conditions, paste0(newnames, "_test_r2"))]
analyzedat_test[, c(1:lc):=lapply(.SD, factor), .SDcols=c(1:lc)]

setwd("./simulatie_2021")
#sink("output.txt") #makes connection with a .txt file, use it later
#unlink('output.txt') #unlinks the connection and deletes .txt file.

analyzedat_test <- na.omit(analyzedat_test)

#Check for significant moderators using eta squared
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

#the dependent variables (the test r2)
yvars <- paste0(newnames, "_test_r2")
#creates a list for every Anova with the results for the anova and the effect sizes for the conditions on the algorithms
anovas<-lapply(yvars, function(yvar){
  form<-paste(yvar, '~(', paste(unlist(conditions[-lc]), "+", collapse = ' '), conditions[lc], ") ^ 2") #the ^2 signifies that we want all possible effects
  thisaov<-aov(as.formula(form), data=analyzedat_test) #change data according to dataframe you are using
  thisetasq<-EtaSq(thisaov)[ , 2]
  list(thisaov, thisetasq)
})

#creates dataframe with effect sizes for all conditions for all algorithms
etasqs<-data.frame(sapply(anovas, function(x){
  tempnums<-x[[2]]
  formatC(tempnums, 2, format="f")
}))

coef.table<-etasqs

#I think these are not necessary
# table.names<-row.names(coef.table)
# #these couple of lines change names of rows to shorter names
# table.names<-gsub("_studies", "", table.names) #this one does nothing because we dont have a variable with '_studies'
# table.names<-gsub("moderators", "M", table.names) #this changes 'moderators' to M
# table.names<-gsub("residual_heterogeneity", "$\\tau^2$", table.names) #nothing
# table.names<-gsub("es", "$\\beta$", table.names) #changes 'es'to '$beta$ but I do not know why
# table.names<-gsub("mean_study_n", "$\\mean{n}$", table.names) #nothing
# table.names<-gsub("model", "Model ", table.names) #capitalizes 'model'

# row.names(coef.table)<-table.names
names(coef.table)<-c("Horseshoe", "Lasso", "Metaforest" ,"RMA")

coeftable.testr2 <- coef.table
#sink("output.txt", append = TRUE)
#sink()

write.csv(coef.table, file =  "EtaSq_testR2.csv")
#sink("output.txt", append = TRUE)
#sink()

#this makes a nice table in .tex file for the coefficients
library(xtable)
print(xtable(coef.table), file="EtaSq_testR2.tex",sanitize.text.function=function(x){x})


#save.image("Analysis_Data.RData")
###### TOT HIER GECHEKT ELI 2021 #####
load("./simulatie_2021/Analysis_Data.RData")


##Figure out which predictors are most important per variable
#function that sorts conditions that have highest eta to lowest
important_predictors <- function(df, column){
  df.sort <- df[order(df[,column], decreasing = T), , drop = FALSE] #sort by particular column
  rnames <- rownames(df.sort) #store row names in variable
  vals <- sort(df.sort[,column], decreasing = TRUE) #store values in variable
  condition <- c() #memory for condition names
  value <- c() #memory for eta square values
  for(i in 1:nrow(df.sort)){
    condition[i] <- rnames[i]
    value[i] <- vals[i]
  }
  return(cbind(condition, value)) #return matrix of condition name with corresponding EtaSq from highest to lowest
}

#for loop to do it for all models
imp_pred_per_model <- list()
for(i in 1:ncol(coef.table)){
  imp_pred_per_model[[i]] <- important_predictors(coef.table, names(coef.table[i]))
}
names(imp_pred_per_model) <- c('Horseshoe', 'Lasso', 'Metaforest', 'RMA')
imp_pred_per_model
#write.csv(imp_pred_per_model, file = "important_predictors_per_model.csv")


######################################
#test R2 per condition per algorithm #
######################################
measure.vars <- names(analyzedat_test)[-c(1:lc)] #test_r2 values
grouping.vars <- quote(list(k_train,
                            mean_n,
                            es,
                            tau2,
                            alpha_mod,
                            moderators,
                            model))

#obtain mean R2 per condition (so mean over 100 datasets)
testR2_per_condition <- analyzedat_test[,lapply(.SD, mean),by=eval(grouping.vars), .SDcols=measure.vars]

#below I want to obtain per condition a rankorder in which algorithm has highest test R2
mean_R2s <- testR2_per_condition[,(lc+1):ncol(testR2_per_condition)] #obtain dataframe without condition names

#obtain maximum value for R2 comparing the algoritms
for(i in 1:nrow(mean_R2s)){
  mean_R2s$highest[i] <- max(mean_R2s[i,1],mean_R2s[i,2],mean_R2s[i,3],mean_R2s[i,4])
}

#if R2 value == highest for that row, then add name to variable algname
for(i in 1:nrow(mean_R2s)){
  if(mean_R2s[i,1] == mean_R2s[i,5]){
    mean_R2s$algnames[i] <- 'Horseshoe'
  } else if(mean_R2s[i,2] == mean_R2s[i,5]){
    mean_R2s$algnames[i] <- 'Lasso'
    } else if(mean_R2s[i,3] == mean_R2s[i,5]){
      mean_R2s$algnames[i] <- 'Metaforest'
    } else if(mean_R2s[i,4] == mean_R2s[i,5]){
      mean_R2s$algnames[i] <- 'RMA'}
  else{stop(paste0('there is no highest value, inspect case ', i))}
}

testR2_per_condition$highestR2 <- mean_R2s$algnames
testR2_per_condition$R2value <- mean_R2s$highest

best_pred_per_condition <- within(testR2_per_condition, rm('hs_test_r2','lasso_test_r2', 'mf_r_test_r2', 'rma_test_r2'))
#write.csv(best_pred_per_condition, file = 'best_performing_algorithm_per_condition.csv')

#save.image("Analysis_Data.RData")
###### TOT HIER GECHEKT ELI 2021 #####
load("./simulatie_2021/Analysis_Data.RData")


##################
# Power analysis #
##################
#subset for test_r2 values
analyzedat <- dat[ , .SD, .SDcols=c(conditions, paste0(newnames, "_test_r2"))]
analyzedatt <- dat[ , .SD, .SDcols=c(conditions, 'mf_r_test_r2', "lasso_test_r2")]

#make conditions factors
analyzedat[, c(1:lc):=lapply(.SD, factor), .SDcols=c(1:lc)]


#for every condition seperately; if the 20% quantile from the test_r2 > 0, then give 1, otherwise
power2<-analyzedatt[,lapply(.SD, function(x){ifelse(quantile(x, probs = .2, na.rm = TRUE)>0, 1, 0)}),by=eval(grouping.vars), .SDcols=measure.vars]
plotdat<-power2

bin_cond <- function()

plotdat[, power := "Neither"] #create variable with value 'Neither'
plotdat[mf_r_test_r2 == 1 & lasso_test_r2 == 0, power := "Only MetaForest"]
plotdat[mf_r_test_r2 == 0 & lasso_test_r2 == 1, power := "Only lasso"]
plotdat[mf_r_test_r2 == 1 & lasso_test_r2 == 1, power := "Both"]
plotdat[, power := factor(power)]
plotdat$power <- ordered(plotdat$power, levels = c("Neither", "Only MetaForest", "Only lasso", "Both"))
table(plotdat$power)

condnames <- c("k", "n", "es", "tau2", "alpha", "M", "Model")
names(plotdat)[c(1:lc)]<- condnames

categorical<-condnames

plotdat[, (categorical) := lapply(.SD, factor), .SDcols=categorical] #make conditions factors
plotdat$tau2 <- factor(plotdat$tau2, labels = c("tau^2:~.00", "tau^2:~.04", "tau^2:~.1")) #change value names for tau2
plotdat$es <- factor(plotdat$es, labels = c("beta:~0","beta:~0.2", "beta:~0.5", "beta:~0.8")) #for effect size

table(analyzedat$model) #there are 4 models, denote them with a, b, c and d
levels(plotdat$Model)<-c("a", "b", "c", "d")

table(plotdat$power)
library(ggplot2)
thisplot <- ggplot(plotdat, aes(x=n, y=k, fill = power)) +
  geom_raster(hjust = 0, vjust = 0)+
  facet_grid(Model+M ~ es+tau2, labeller = labeller(es = label_parsed, tau2 = label_parsed, Model=label_both, M = label_both))+
  scale_fill_manual(values=c("white", "grey50", "grey25", "black"))+
  scale_x_discrete(expand=c(0,0))+
  scale_y_discrete(expand=c(0,0))+
  labs(x=expression(bar(n)), fill="Power > .80")
#theme(legend.title = element_blank()) +
#thisplot+annotate("text", label = "p[mf<mc]", parse=TRUE)

#have to adapts width and height
ggsave("power_plot.pdf", plot=thisplot, device="pdf", width=210, height=297, units="mm")

#Percentage table
power <- analyzedat[!(model==1),lapply(.SD, function(x){sum(x > 0)/length(x)}),by=eval(grouping.vars), .SDcols=measure.vars]
powertmp <- analyzedat[model==1,lapply(.SD, function(x){sum(x < 0)/length(x)}),by=eval(grouping.vars), .SDcols=measure.vars]
power <- rbindlist(list(powertmp, power))

power <- power
power <- melt(power, measure.vars = c("mf_r_test_r2", "lasso_test_r2")) #melt() makes wide to long format for data.table, but
#only for the measure.vars

tmp <- dcast(power, model+moderators+tau2~variable+ es + k_train + mean_n + alpha_mod) #dcast() is wide to long
#the conditoins on the left side of the formula remain in long format, while the conditions on the right side will
#be put in variable names and thus wide format
write.csv(tmp, "study1 power.csv") #writes tmp to excel file .csv format

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
#table.names<-gsub("_studies", "", table.names)
#table.names<-gsub("moderators", "M", table.names)
#table.names<-gsub("residual_heterogeneity", "$\\tau^2$", table.names)
#table.names<-gsub("es", "$\\beta$", table.names)
#table.names<-gsub("mean_study_n", "$\\mean{n}$", table.names)
#table.names<-gsub("model", "Model ", table.names)

row.names(coef.table)<-table.names
names(coef.table)<-c("Horseshoe", "Lasso", "Metaforest", "RMA")


write.csv(coef.table, file =  "EtaSq_testR2.csv")
#sink("output.txt", append = TRUE)
#sink()

library(xtable)
print(xtable(coef.table), file="EtaSq_testR2.tex",sanitize.text.function=function(x){x})


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
