#save.image('./simulatie_2021/testfile.RData')
#load('./simulatie_2021/testfile.RData')

library(readxl)
library(metafor)
library(metaforest)
library(pema)
library(rstan)
library(brms)
#library(rstantools)
library(bayesplot)


#load in data that caspar provided and get a sample of 100 rows
data2 <- read.csv('./simulatie_2021/Data_SJ.csv') #read from csv
dat.dup <- data2[,-1] #exclude first variable as that is index indicator
dat <- unique(dat.dup)
rm(data2, dat.dup)

#lets reduce the dataset down a bit to contain only the study ID's that are duplicates
datdups <- rle(dat$study_ID)
dupsin <- which(datdups$lengths>1)
datsub_mes <- dat[dupsin,]
rm(datdups, dupsin)

#try converting grouping factor study to integer
#make it a little smaller
set.seed(6164900)
datsub_mes2 <- datsub_mes[sample(nrow(datsub_mes), 100), ]
library(dplyr)
datsub_mes2 <- datsub_mes2 %>% mutate(group_no = as.integer(factor(study_ID)))
datsub_mes2[2,c('discussion', 'published', 'dilemma_type')] <- NA #to check if everything works with NA values


#run rma
rma.re <- rma(yi, vi, mods = cbind(discussion, dilemma_type), data = datsub_mes2) #works
summary(rma.re)

#run random forest
datsub_mes2 <- missRanger::missRanger(datsub_mes2)
fit.mf <- MetaForest(yi ~ discussion + published + dilemma_type, num.trees = 150, data = datsub_mes2) #does not work with NA values
summary(fit.mf)

#fit brma models #also does not work with missings
fit.lasso <- pema::brma(yi ~ discussion + dilemma_type,
                        study = datsub_mes2$group_no,
                        vi = datsub_mes2$vi,
                        data = datsub_mes2,
                        method = 'lasso')
summary(fit.lasso)
#With text vi
fit.lasso <- pema::brma(yi ~ discussion + dilemma_type,
                        study = datsub_mes2$group_no,
                        vi = "vi",
                        data = datsub_mes2,
                        method = 'lasso')

fit.hs <- pema::brma(yi ~ discussion + dilemma_type,
                        study = datsub_mes2$group_no,
                        vi = datsub_mes2$vi,
                        data = datsub_mes2,
                        method = 'hs')
summary(fit.hs)








#
fit.lasso2 <- pema::brma(yi ~ discussion + dilemma_type,
                        study = as.numeric(as.factor(datsub_mes2$group_no)),
                        vi = datsub_mes2$vi,
                        data = datsub_mes2,
                        method = 'lasso')

