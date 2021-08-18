
###################################
# Predictive performance analyses #
###################################
library(data.table)
dat <- as.data.table(readRDS(r"(C:\Users\e_lib\OneDrive\Documents\Caspar Repositories\simulatie_2021\sim_results_2021-08-11.RData)"))
newnames <- c('hs', 'lasso', 'mf_r', 'rma') #names of models
conditions <- c("k_train", "mean_n", "es", "tau2","alpha_mod", "moderators", "model") #conditions that vary
lc <- length(conditions) #length of conditions, handy for further code

#subsets the data based on particular performance criterion
#in this case train_r2
analyzedat <- dat[ , .SD, .SDcols=c(conditions, paste0(newnames, "_train_r2"))]

#converts the condition variables to factors
analyzedat[, c(1:lc):=lapply(.SD, factor), .SDcols=c(1:7)]

setwd("./simulatie_2021")
#sink("output.txt") #makes connection with a .txt file, use it later
#unlink('output.txt') #unlinks the connection and deletes .txt file.


#there are missing values in the data
apply(analyzedat, 2, function(x){sum(is.na(x))}) #17 missing values in rma_train_r2
cases_rma_mis <- which(is.na(analyzedat$rma_train_r2)) #which cases are missing
rma_mis <- analyzedat[cases_rma_mis,] #subset data for missing cases
summary(rma_mis) #at first sight, I do not see an obvious pattern to why these cases are missing

#Since only 17 cases out of 4000000 are missing, I do not think it will have a big impact on further analyses
#Although it might be interesting to see why the cases are missing, for now I will continue with the observed data
analyzedat <- na.omit(analyzedat)

#obtain mean and standard deviations
tmp<-rbind(round(analyzedat[model != "1", lapply(.SD, mean), .SDcols=c((lc+1):ncol(analyzedat))], 2),
           round(analyzedat[model !="1", lapply(.SD, sd), .SDcols=c((lc+1):ncol(analyzedat))], 2))
tmp<-data.frame(names(tmp), t(tmp))
apply(tmp, 1, function(x){paste(x, collapse=", ")})
colnames(tmp) <- c('alg', 'mean', 'sd')
View(tmp) #r2 does not seems to make much sense for hs and lasso
sink()

#this functions plots densities for the chains and provides descriptives
pre_analysis <- function(subdat = analyzedat, lb = -Inf, ub = Inf){ #lb = lowerbound, ub = upperbound for the boundaries of impossible values
  #plot densities for chains in one window
  par(mfrow = c(2, 2))
  subdat <- as.data.frame(subdat)
  for(column in (lc+1):ncol(subdat)){
    col_obs <- na.omit(subdat[,column]) #there are sometimes missings so omit them
    prop.impvals <- sum(col_obs < lb | col_obs > ub) / length(col_obs) #proportion of impossible values
    plot(density(col_obs),
         main = paste0('histogram of ', colnames(subdat)[column]),
         xlab = paste0('proportion impossible values: ', round(prop.impvals, 3)))
  }
  print(psych::describe(subdat[,(lc+1):ncol(subdat)])) #print descriptives of chains
}

pre_analysis(analyzedat, lb = 0 , ub = 1) #chains fo rma and mf look fine, not so much for hs and lasso.
#also, both hs and lasso have ~60% of r2 value that are impossible to obtain.




########################################
#Predictive performance on testing data#
########################################

analyzedat <- dat[ , .SD, .SDcols=c(conditions, paste0(newnames, "_test_r2"))]

#converts the condition variables to factors
analyzedat[, c(1:lc):=lapply(.SD, factor), .SDcols=c(1:lc)]

setwd("./simulatie_2021")
#sink("output.txt") #makes connection with a .txt file, use it later
#unlink('output.txt') #unlinks the connection and deletes .txt file.


#there are missing values in the data
apply(analyzedat, 2, function(x){sum(is.na(x))}) #17 missing values in rma_test_r2 as well
cases_rma_mis2 <- which(is.na(analyzedat$rma_test_r2)) #which cases are missing
sum(cases_rma_mis != cases_rma_mis2) #exacly the same cases are missing for rma_train_r2 as for rma_test_r2

analyzedat <- na.omit(analyzedat)

#obtain mean and standard deviations
tmp<-rbind(round(analyzedat[model != "1", lapply(.SD, mean), .SDcols=c((lc+1):ncol(analyzedat))], 2),
           round(analyzedat[model !="1", lapply(.SD, sd), .SDcols=c((lc+1):ncol(analyzedat))], 2))
tmp<-data.frame(names(tmp), t(tmp))
apply(tmp, 1, function(x){paste(x, collapse=", ")})
colnames(tmp) <- c('alg', 'mean', 'sd')
View(tmp) #r2 does not seems to make much sense for hs and lasso
sink()

pre_analysis(analyzedat, lb = 0 , ub = 1) #same problem,
#both hs and lasso have ~55% of r2 value that should be impossible to obtain.
#although both mf and rma also contain impossible values.

###checked until here, Eli 2021
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

yvars<-c(paste0(newnames, "_test_r2"), paste(newnames[1], newnames[-1], sep = "_")) #I do not really get how this helps
yvars <- paste0(newnames, "_test_r2")


withoutmodel1<-analyzedat[model!="1",] # I do not think this is necessary
fittest <- aov(hs_test_r2 ~ (k_train + es) ^ 2, analyzedat)#the ^2 signifies that we want all possible effects
summary(fittest)
EtaSq(fittest)

##lets first do the anovas for a subsample of the data to check how and if everything works
set.seed(6164900)
subdat <- analyzedat[sample(nrow(analyzedat), 100), ]

#creates a list for every Anova with the results for the anova and the effect sizes for the conditions on the algorithms
anovas<-lapply(yvars, function(yvar){
  form<-paste(yvar, '~(', paste(unlist(conditions[-lc]), "+", collapse = ' '), conditions[lc], ") ^ 2")
  thisaov<-aov(as.formula(form), data=subdat) #change data according to dataframe you are using
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
#study1.coef.table.withoutmodel1
#sink()

#study1.coef.table.withoutmodel1 <- study1.coef.table.withoutmodel1[, -c(2,3, 7, 8)]

#this makes a nice table in .tex file for the coefficients
library(xtable)
print(xtable(coeftable.testr2), file="coeftable.testr2.tex",sanitize.text.function=function(x){x})


##Figure out which predictors are most important per variable
as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}
#sink("output.txt", append = TRUE)
lapply(coeftable.testr2, function(x){
  tmp <- as.numeric.factor(x) #as.numeric.factor() is not an existing function
  names(tmp)<- names(coeftable.testr2)
  sort(tmp, decreasing = TRUE)
})
#sink()

#functoin that sorts conditions that have highest eta to lowest
important_predictors <- function(df, column){
  df.sort <- df[order(df[,column]), , drop = FALSE]
  rnames <- rownames(df.sort)
  vals <- sort(df.sort[,column], decreasing = TRUE)
  memory <- c()
  for(i in 1:nrow(df.sort)){
    elem <- paste0(i, ".) ", rnames[i], "= ", vals[i])
    memory[i] <- elem
  }
  return(memory)
}

#for loop to do it for all models
imp_pred_per_model <- list()
for(i in 1:ncol(coeftable.testr2)){
  imp_pred_per_model[[i]] <- important_predictors(coeftable.testr2, names(coeftable.testr2[i]))
}
names(imp_pred_per_model) <- c('Horseshoe', 'Lasso', 'Metaforest', 'RMA')
imp_pred_per_model

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
