# 8/27/20
# Exponential model analysis for Example \ref{ex:ExpModelAnalysis}=11.3
library("survival")
ds <- read.csv("CALGB8541.csv",header=TRUE)

dsarm1 <- ds[ds$arm==1,]
n <- nrow(dsarm1)
nu <- n-sum(dsarm1$dfsstat)
nybar <- sum(dsarm1$dfsyrs)

a0 <- 1.4
b0 <- 5
lambda <- rgamma(10000,a0,b0)
Einvlambda <- mean(1/lambda) 
ll <- 1/qgamma(0.975,a0,b0)
uu <- 1/qgamma(0.025,a0,b0)
Einvlambda # prior expectation of mean disease-free survival
ll         # prior lower limit for mean disease-free survival (95% probI)
uu         # prior lower limit for mean disease-free survival (95% probI)
# Results ############
## > [1] 11.50825
## > [1] 1.114608
## > [1] 57.59722
######################

a1 <- a0+nu
b1 <- b0+nybar
postlambda <- rgamma(10000,a1,b1)
Einvpostlambda <- mean(1/postlambda)
postinterval95 <- quantile(1/postlambda,probs=c(0.025,0.975))
Einvpostlambda
postinterval95
# Results ############
## > [1] 16.70207     # posterior expected mean disease-free survival
## >  
##   2.5%    97.5% 
## 14.91289 18.67308  # posterior 95% probI for mean disease-free survival
######################

pdf("ExpModelSurvSamples.pdf")
t <- seq(0.1,20,0.1)
plot(t,exp(-postlambda[1]*t),ylim=c(0,1),ylab="Posterior Samples of S(t)", type="l",lwd=0.5)
for (i in 2:10) {lines(t,exp(-postlambda[i]*t),type="l",lwd=0.5)}
dev.off()

pdf("ExpModelAnalysis.pdf")
dfsobject <- Surv(dsarm1$dfsyrs,1-dsarm1$dfsstat)
plot(dfsobject,conf.int=FALSE,ylab="Disease-free survival function",xlab="Disease-free survival")
meanS <- rep(0,length(t))
for(i in 1:10000){meanS <- meanS+exp(-postlambda[i]*t)}
meanS <- meanS/10000
lines(t,meanS,type="l")
lines(t,exp(-as.numeric(1/postinterval95[1])*t),type="l",lty=2)
lines(t,exp(-as.numeric(1/postinterval95[2])*t),type="l",lty=2)
dev.off()
