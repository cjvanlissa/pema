# Based on: http://127.0.0.1:21883/library/glmmLasso/demo/glmmLasso-soccer.r
################## More Elegant Method ############################################
## Idea: start with big lambda and use the estimates of the previous fit (BUT: before
## the final re-estimation Fisher scoring is performed!) as starting values for the next fit;
## make sure, that your lambda sequence starts at a value big enough such that all covariates are
## shrinked to zero;

## Using BIC (or AIC, respectively) to determine the optimal tuning parameter lambda
lambda <- seq(2,0, length.out = 20)

glmmLasso(fix=yi ~ 1, rnd = ~ 1 | study, data = dat, lambda)
dat$studfac <- as.factor(dat$study)
X <- model.m
res.glml <- glmmLasso(fix=yi ~ 1 + as.factor(controls) + as.factor(design) + as.factor(a_measure)+ as.factor(c_measure) + meanage + quality, rnd = list(studfac = ~ 1), data = dat, lambda = 1)

BIC_vec<-rep(Inf,length(lambda))

# specify starting values for the very first fit; pay attention that Delta.start has suitable length!
Delta.start<-as.matrix(t(rep(0,dim(res.glml$Deltamatrix)[2])))
Q.start<-0.1

for(j in 1:length(lambda)){
  print(paste("Iteration ", j,sep=""))

  res <- glmmLasso(fix=yi ~ 1 +
                     as.factor(controls) +
                     as.factor(design) +
                     as.factor(a_measure)+
                     as.factor(c_measure) +
                     meanage +
                     quality,
                   rnd = list(studfac = ~ 1),
                   data = dat,
                   lambda=lambda[j],
                   switch.NR=F,
                   final.re=TRUE,
                   control = list(start=Delta.start[j,],q_start=Q.start[j]))

  print(colnames(res$Deltamatrix)[2:7][res$Deltamatrix[res$conv.step,2:7]!=0])
  BIC_vec[j]<-res$bic
  Delta.start<-rbind(Delta.start,res$Deltamatrix[res$conv.step,])
  Q.start<-c(Q.start,res$Q_long[[res$conv.step+1]])
}

opt3<-which.min(BIC_vec)

res_final <- glmmLasso(fix=yi ~ 1 +
                          as.factor(controls) +
                          as.factor(design) +
                          as.factor(a_measure)+
                          as.factor(c_measure) +
                          meanage +
                          quality,
                        rnd = list(studfac = ~ 1),
                        data = dat,
                       lambda=lambda[opt3],
                       switch.NR=F,final.re=TRUE,
                       control = list(start=Delta.start[opt3,],q_start=Q.start[opt3]))


summary(res_final)

## plot coefficient paths
par(mar=c(6,6,4,4))
plot(lambda,Delta.start[2:(length(lambda)+1),2],type="l",ylim=c(-1e-1,1e-1),ylab=expression(hat(beta[j])))
lines(c(-1000,1000),c(0,0),lty=2)
for(i in 3:7){
  lines(lambda[1:length(lambda)],Delta.start[2:(length(lambda)+1),i])
}
abline(v=lambda[opt3],lty=2)

