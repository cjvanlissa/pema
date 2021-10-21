#creates traceplots for metrics
Traceplot <- function(test, train, param, alg){
  plot(x = 1:length(test), y = test, type = 'n', main = paste('Traceplot Test and Train', param , ' for', alg), xlab = 'iterations', ylab = paste('values for ', param)) #empty plot with correct dimensions for the traced values
  lines(x = 1:length(test),  y = test, col = 'red') #superimpose traced values for testing r2
  lines(x= 1: length(train), y = train, col = 'blue') #superimpose traced values for training r2
  legend('bottomleft', c('Train R2', 'Test R2'), lty = 1, col = c('blue', 'red'), bg = 'white', cex = 0.75) #add informative legend
}


#creates densities for metric
plotdens <- function(df){
  df <- as.data.frame(df)
  for(column in (lc+1):ncol(df)){ #plot densities for all algorithms
    plot(density(df[,column]),
         main = paste0('histogram of ', colnames(df)[column])
         #xlim = c(-30,1)
    )
  }
}

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

#functions for plotting marginal and conditional effects of design factors for algorithms
plot_marginal_mean <- function(condition, df, lonames, metric){
  condition <- as.name(condition)
  hs <- as.name(lonames[1])
  las <- as.name(lonames[2])
  mf <- as.name(lonames[3])
  rma <- as.name(lonames[4])

  cond <- df %>%
    group_by(!!condition)%>%
    summarise(Horseshoe = mean(!!hs),
              Lasso = mean(!!las),
              MetaForest = mean(!!mf),
              RMA = mean(!!rma))%>%
    pivot_longer(cols = !all_of(condition), names_to = 'alg', values_to = metric)

  ggplot(data = cond, aes(x = !!condition,y = !!as.name(metric), linetype = alg, group = alg, shape = alg)) +
    geom_line(size = 1, ) +
    geom_point(size = 2.5) +
    theme_classic(base_size = 25)
}

plot_interaction_mean <- function(condition, intcond, df, lonames, metric){
  condition <- as.name(condition)
  intcond <- as.name(intcond)
  hs <- as.name(lonames[1])
  las <- as.name(lonames[2])
  mf <- as.name(lonames[3])
  rma <- as.name(lonames[4])

  cond <- df %>%
    group_by(!!condition, !!intcond)%>%
    summarise(Horseshoe = mean(!!hs),
              Lasso = mean(!!las),
              Metaforest = mean(!!mf),
              RMA = mean(!!rma)) %>%
    pivot_longer(cols = c(Horseshoe, Lasso, Metaforest, RMA), names_to = 'alg', values_to = metric)

  ggplot(data = cond, aes(x = !!condition,y = !!as.name(metric), linetype = alg, group = alg, shape = alg)) +
    geom_line(size = 1, ) +
    geom_point(size = 2.5) +
    facet_wrap(~cond[[intcond]], scales = 'free') +
    theme_classic(base_size = 25)
}

plot_marginal_median <- function(condition, df, lonames, metric){
  condition <- as.name(condition)
  hs <- as.name(lonames[1])
  las <- as.name(lonames[2])
  mf <- as.name(lonames[3])
  rma <- as.name(lonames[4])

  cond <- df %>%
    group_by(!!condition)%>%
    summarise(Horseshoe = median(!!hs),
              Lasso = median(!!las),
              Metaforest = median(!!mf),
              RMA = median(!!rma)) %>%
    pivot_longer(cols = !all_of(condition), names_to = 'alg', values_to = metric)

  ggplot(data = cond, aes(x = !!condition,y = !!as.name(metric), linetype = alg, group = alg, shape = alg)) +
    geom_line(size = 1, ) +
    geom_point(size = 2.5) +
    theme_classic(base_size = 25)
}

plot_interaction_median <- function(condition, intcond, df, lonames, metric){
  condition <- as.name(condition)
  intcond <- as.name(intcond)
  hs <- as.name(lonames[1])
  las <- as.name(lonames[2])
  mf <- as.name(lonames[3])
  rma <- as.name(lonames[4])

  cond <- df %>%
    group_by(!!condition, !!intcond)%>%
    summarise(Horseshoe = median(!!hs),
              Lasso = median(!!las),
              Metaforest = median(!!mf),
              RMA = median(!!rma)) %>%
    pivot_longer(cols = c(Horseshoe, Lasso, Metaforest, RMA), names_to = 'alg', values_to = metric)

  ggplot(data = cond, aes(x = !!condition,y = !!as.name(metric), linetype = alg, group = alg, shape = alg)) +
    geom_line(size = 1, ) +
    geom_point(size = 2.5) +
    facet_wrap(~cond[[intcond]], scales = 'free') +
    theme_classic(base_size = 25)
}



#this functions creates a subsetted dataframe for the wanted sel column and the number of desired moderators
# and calculates true and false positives an negatives
library(tidyverse)
sel_subset <- function(moderators, pattern){
  conditions2 <- c(conditions, 'identifier')
  number <- as.numeric(moderators)
  varsel_sub <- varsel[varsel$moderators == number,] #subset for number of desired moderators

  #subset for desired columns for example hs_ci or rma_ci
  sel_sub <- str_extract(colnames(varsel_sub), pattern = pattern)
  sel_sub <-sel_sub[!is.na(sel_sub)]
  varsel_sub <- varsel_sub[ , .SD, .SDcols=c(conditions2, sel_sub)]

  #finally omit columns containing only NA's
  non_empty_cols <- (colSums(is.na(varsel_sub) | varsel_sub == "") != nrow(varsel_sub)) #extract non-empty columns
  varsel_sub <- tidyr::as_tibble(varsel_sub)
  varsel_sub <- varsel_sub[ ,non_empty_cols]

  #subset to null model and es != 0
  null_model <- varsel_sub[varsel_sub$es == 0, ]
  sig_model <- varsel_sub[varsel_sub$es != 0, ]

  #obtain fractions true negatives, false positives for null_model
  null_model$TN <- apply(null_model[, (lc+2):(ncol(null_model))], 1, function(x){sum(x == 0)/moderators})
  #null_model$FP <- 1- null_model$TN

  #obtain for sig model and also true positives and false negatives
  sig_model$TP <- apply(sig_model[, (lc+2)], 1, function(x){sum(x == 1)})
  sig_model$TN <- apply(sig_model[, (lc+3):(ncol(sig_model) - 1)], 1, function(x){sum(x==0)/(moderators-1)})
  #sig_model$FP <- 1- sig_model$TN
  #sig_model$FN <- 1- sig_model$TP
  #sig_model$FPS <- FPSf(sig_model$TP, sig_model$TN)

  #join the two models together
  full_model <- suppressMessages(full_join(sig_model, null_model))
  full_model <- subset(full_model, select=-c(identifier))
  full_model <- as.data.table(full_model)

  return(full_model[])
}



#this functions binds the dfs per algorithm for all levels of moderators
sel_all_levels <- function(pattern){

  # extract names which are used later
  metric_name <- str_extract(pattern, pattern = "[:alpha:]+_[:alpha:]+")
  algname <- str_extract(metric_name, pattern = "[:alpha:]+")

  #create memory and run sel_subset function for all number of moderators
  df <- c()
  for(i in as.numeric(names(table(varsel$moderators)))){
    metrics <- sel_subset(i, pattern)
    df <- vctrs::vec_c(df, metrics)
    rm(metrics) }

  #give appropriate names and subset for necessary column
  names(df)[grep("(TP|TN)", colnames(df))] <- paste(metric_name, names(df)[grep("(TP|TN)", colnames(df))], sep = "_")
  irr_colnums <- grep(paste0("^",algname,".*\\d$"), colnames(df))
  irr_colnames <- colnames(df)[irr_colnums]
  df <- df[, -..irr_colnames]
  df$identifier <- 1:nrow(df)
  return(df)

}
