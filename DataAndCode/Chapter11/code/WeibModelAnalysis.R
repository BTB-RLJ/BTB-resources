# 9/8/20
# Posterior analysis and plots follow prior analysis
# gibbsWeibull() tested with default simulated data; simulation fixed!

# Weibull model analysis for Example \ref{ex:WeibModelAnalysis}=11.4
library("survival")
ds <- read.csv("CALGB8541.csv",header=TRUE)

d <- 3
c <- 1.5  # follows from P(alpha<1)=0.5
alphastar <- 1.125  # E(alpha)=1.125
tmedian <- 12   # educated guess for median survival time
lambdastar <- log(2)/tmedian^alphastar
u95 <- 20
lambda05 <- log(2)/u95^alphastar # 5th percentile of lambda
# trial and error to find b such that lambda ~ Ga(1+lambdastar*b,b) has
# 5th percentile lambda05.
b <- 3
qgamma(0.05,1+lambdastar*b,b)-lambda05

## > b <- 3
## > qgamma(0.05,1+lambdastar*b,b)-lambda05
## [1] 0.001674735
## > b <- 4
## > qgamma(0.05,1+lambdastar*b,b)-lambda05
## [1] -0.002336409
## > b <- 3.5
## > qgamma(0.05,1+lambdastar*b,b)-lambda05
## [1] -0.0006355157
## > b <- 3.25
## > qgamma(0.05,1+lambdastar*b,b)-lambda05
## [1] 0.0004258248
## > b <- 3.375
## > qgamma(0.05,1+lambdastar*b,b)-lambda05
## [1] -0.000125682
## > b <- 3.325
## > qgamma(0.05,1+lambdastar*b,b)-lambda05
## [1] 8.965625e-05
## >

a <- 1+lambdastar*b

## > a <- 1+lambdastar*b
## > a
## [1] 1.140779
## > b
## [1] 3.325


# Prior simulation; same for both outcomes
a <- 1.141
b <- 3.325
c <- 1.5
d <- 3
phi <- c/runif(5000)^(1/d)
alpha <- runif(5000,0,phi)
lambda <- rgamma(5000,a,b)
mediansurv <- log(2)/lambda^(1/alpha)

# Prior means and intervals

alphaPriorE <- mean(alpha)
alphaPrior95 <- quantile(alpha,probs=c(0.025,0.975))
lambdaPriorE <- mean(lambda)
lambdaPrior95 <- quantile(lambda,probs=c(0.025,0.975))
mediansurvPriorMedian <- median(mediansurv)
mediansurvPrior95 <- quantile(mediansurv,probs=c(0.025,0.975))

alphaPriorE
## > alphaPriorE
## [1] 1.130105

alphaPrior95
## > alphaPrior95
## > alphaPrior95
##       2.5%      97.5% 
## 0.04666052 3.38354939 

lambdaPriorE
## > lambdaPriorE
## [1] 0.344953

lambdaPrior95
## > lambdaPrior95
##       2.5%      97.5% 
## 0.01317663 1.18118727

mediansurvPriorMedian
## > mediansurvPriorMedian
## [1] 2.889837

mediansurvPrior95
## > mediansurvPrior95
##         2.5%        97.5% 
## 5.708669e-01 4.102470e+12

# Prior information translated to survival functions
t <- seq(0.1,20,0.1)     # time grid for plotting
survatgrid <- matrix(rep(0,length(lambda)*length(t)),length(lambda),length(t))
for(i in 1:length(lambda)){survatgrid[i,] <- exp(-lambda[i]*t^alpha[i])}
survPriorMedian <- apply(survatgrid,2,"median")
survPrior95 <- apply(survatgrid,2,"quantile",probs=c(0.025,0.975))

# Plot prior mean survival function and intervals
plot(t,survPriorMedian,type="l")
lines(t,survPrior95[1,],type="l",lty=2)
lines(t,survPrior95[2,],type="l",lty=2)
    

# Posterior sampling added 9/8/20 

dsarm1 <- ds[ds$arm==1,]

# code from Chapter 4, modified to allow rt-censoring
# The changes are: delta as input and logprody now
# over delta=1; n replaced by n_u.


require("HI")
logdensity <- function(alpha,phi,y,n_u,logprody,lambda)
  (n_u*log(alpha)+alpha*logprody-lambda*sum(exp(alpha*log(y))))
# censoring change: n to n_u
gibbsWeibull <- function(
  y=NULL,                # default data from Weib(alpha=2,lambda=1)
  delta=NULL,            # default rt-censoring with Exp(1)
  a=1.0,                 # shape parameter in prior for lambda
  b=1.0,                 # rate parameter in prior for lambda
  c=1,                   # Pareto lower bound parameter
  d=2,                   # Pareto shape parameter
  burnin=100,            # initial iterations to discard
  thin=1,                # thining factor; no thining by default
  M=1000,                # number of samples returned
  lambdainit=NULL,       # initial value for lambda
  alphainit=NULL,        # initial value for Weibull shape parameter
  phiinit=NULL           # initial value for hierarchical parameter
  )
{
  if(length(y)==0) {t <- rexp(500)^0.5
                    hammer <- rexp(500)
                    y <- apply(cbind(t,hammer), 1, "min")
                    delta <- (t<hammer)
  }
  if(length(lambdainit)==0) lambdainit <- a/b
  if(length(phiinit)==0) phiinit <- c*d/(d-1)
  if(length(alphainit)==0) alphainit <- phiinit/2
  lambdas <- rep(0,M);  alphas <- rep(1,M)
  lambda <- lambdainit; alpha <- alphainit; phi <- phiinit
  n_u <- sum(delta);  logprody <- sum(log(y)*delta) # censoring change
  a1 <- a+n_u                                       # censoring change
  d1 <- d+1
  for(g in (1-burnin):(M*thin)){ # begin Gibbs loop
    b1 <- b+sum(exp(alpha*log(y)))
    lambda <- rgamma(1,shape=a1,rate=b1)
    phi <- max(c,alpha)/runif(1)^(1/d1)
    alpha <- arms(alpha,logdensity,
                  function(x,...) ((x>0)*(x<phi)),1,phi,
                  y,n_u,logprody,lambda)            # censoring change
    if(g>0 & (g%%thin)==0){      # save posterior sample
      lambdas[g] <- lambda;   alphas[g] <- alpha
    }
  }                              # end Gibbs loop
  return(list(lambda=lambdas,alpha=alphas))
}


############################
# overall-survival outcome
############################

survsamples <- gibbsWeibull(dsarm1$survyrs,1-dsarm1$survstat,
                            a=a,b=b,c=c,d=d)

summary(survsamples$lambda)
## > summary(survsamples$lambda)
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
## 0.01724 0.02763 0.03161 0.03174 0.03570 0.05188 

summary(survsamples$alpha)
## > summary(survsamples$alpha)
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##  0.9426  1.0783  1.1230  1.1265  1.1681  1.3452 

mediansurv <- log(2)/survsamples$lambda^(1/survsamples$alpha)

# Posterior means and intervals

alphaPostE <- mean(survsamples$alpha)
alphaPost95 <- quantile(survsamples$alpha,probs=c(0.025,0.975))
lambdaPostE <- mean(survsamples$lambda)
lambdaPost95 <- quantile(survsamples$lambda,probs=c(0.025,0.975))
mediansurvPostE <- mean(mediansurv)
mediansurvPost95 <- quantile(mediansurv,probs=c(0.025,0.975))

alphaPostE
## > alphaPost95
##     2.5%    97.5% 
## 1.005169 1.260631 

alphaPost95
## > alphaPost95
##     2.5%    97.5% 
## 1.005169 1.260631 

lambdaPostE
## > lambdaPostE
## [1] 0.031744

lambdaPost95
## > lambdaPost95
##       2.5%      97.5% 
## 0.02138283 0.04345342 

mediansurvPostE
## > mediansurvPostE
## [1] 15.1009

mediansurvPost95
## > mediansurvPost95
##     2.5%    97.5% 
## 13.45921 17.19884 


# Posterior information translated to survival functions
t <- seq(0.1,20,0.1)     # time grid for plotting
survatgrid <- matrix(rep(0,length(survsamples$lambda)*length(t)),length(survsamples$lambda),length(t))
for(i in 1:nrow(survatgrid)){survatgrid[i,] <- exp(-survsamples$lambda[i]*t^survsamples$alpha[i])}
survPostE <- apply(survatgrid,2,"mean")
survPost95 <- apply(survatgrid,2,"quantile",probs=c(0.025,0.975))
    
# Plot (prior too?)

pdf("WeibModelAnalysisOvSurv.pdf")
Survobject <- Surv(dsarm1$survyrs,1-dsarm1$survstat)
plot(Survobject,conf.int=FALSE,ylab="Overall survival probability",xlab="Overall survival time",lwd=3,lty=3)
lines(t,survPostE,type="l",lwd=3,lty=1)
lines(t,survPost95[1,],type="l",lty=1,lwd=1)
lines(t,survPost95[2,],type="l",lty=1,lwd=1)

lines(t,survPriorMedian,type="l",lty=2,lwd=3)
lines(t,survPrior95[1,],type="l",lty=2,lwd=1)
lines(t,survPrior95[2,],type="l",lty=2,lwd=1)

legend(10,0.9,legend=c("Prior","Posterior","K-M estimate"),lwd=3,lty=c(2,1,3))
dev.off()


###############################
# disease-free-survival outcome
###############################

dfssamples <- gibbsWeibull(dsarm1$dfsyrs,1-dsarm1$dfsstat,a=a,b=b,c=c,d=d)

summary(dfssamples$lambda)
## > summary(dfssamples$lambda)
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
## 0.05574 0.07654 0.08297 0.08319 0.08957 0.11142 

summary(dfssamples$alpha)
## > summary(dfssamples$alpha)
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##  0.7596  0.8422  0.8672  0.8694  0.8950  1.0027 

mediansurv <- log(2)/dfssamples$lambda^(1/dfssamples$alpha)

# Posterior means and intervals

alphaPostE <- mean(dfssamples$alpha)
alphaPost95 <- quantile(dfssamples$alpha,probs=c(0.025,0.975))
lambdaPostE <- mean(dfssamples$lambda)
lambdaPost95 <- quantile(dfssamples$lambda,probs=c(0.025,0.975))
mediansurvPostE <- mean(mediansurv)
mediansurvPost95 <- quantile(mediansurv,probs=c(0.025,0.975))

alphaPostE
## > alphaPostE
## [1] 0.8693628

alphaPost95
## > alphaPost95
##      2.5%     97.5% 
## 0.7959585 0.9534136 

lambdaPostE
## > lambdaPostE
## [1] 0.08318812

lambdaPost95
## > lambdaPost95
##       2.5%      97.5% 
## 0.06540683 0.10266312 

mediansurvPostE
## > mediansurvPostE
## [1] 12.24115

mediansurvPost95
## > mediansurvPost95
##     2.5%    97.5% 
## 10.64397 14.05519 

# Posterior information translated to survival functions
t <- seq(0.1,20,0.1)     # time grid for plotting
survatgrid <- matrix(rep(0,length(dfssamples$lambda)*length(t)),length(dfssamples$lambda),length(t))
for(i in 1:nrow(survatgrid)){survatgrid[i,] <- exp(-dfssamples$lambda[i]*t^dfssamples$alpha[i])}
survPostE <- apply(survatgrid,2,"mean")
survPost95 <- apply(survatgrid,2,"quantile",probs=c(0.025,0.975))
    
pdf("WeibModelAnalysisDfSurv.pdf")
Survobject <- Surv(dsarm1$dfsyrs,1-dsarm1$dfsstat)
plot(Survobject,conf.int=FALSE,ylab="Disease-free survival probability",xlab="Disease-free survival time",lwd=3,lty=3)

lines(t,survPostE,type="l",lwd=3,lty=1)
lines(t,survPost95[1,],type="l",lty=1,lwd=1)
lines(t,survPost95[2,],type="l",lty=1,lwd=1)

lines(t,survPriorMedian,type="l",lty=2,lwd=3)
lines(t,survPrior95[1,],type="l",lty=2,lwd=1)
lines(t,survPrior95[2,],type="l",lty=2,lwd=1)

legend(10,0.9,legend=c("Prior","Posterior","K-M estimate"),lwd=3,lty=c(2,1,3))
dev.off()

