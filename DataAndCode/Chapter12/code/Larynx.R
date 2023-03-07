# remove.packages("DPWeibull")
# setwd("C:/Users/ys8wp/Desktop/")
library(DPWeibull)
library(survival)

LarynxData <- read.csv("LarynxData.csv",header=TRUE)
## str(LarynxData)

## DPresult <- dpweib(Surv(ObsTime,EventIndicator)~factor(Stage)+Age+YrDx,
                   data=LarynxData,
                   predtime = seq(from=0,to=15,length=200))

## newdata <- NULL
## newdata$Stage <-c (1,2,3,4)
## newdata$Age <- rep(65,4)
## newdata$YrDx <- rep(75,4)
## newdata <- data.frame(newdata)
## DPpredict <- predict(DPresult,newdata)

load("DPresult.RData")
load("DPpredict.RData")

#pdf("Larynx.pdf")
# at time 0, log hazard ratio is undefined
pdf("LdataStage4v1LogHR.pdf")
plot(DPresult$predtime,DPresult$loghr.est[3,],
     main="Log hazard ratio for Stage 4 vs Stage 1",
     type="l",lwd=4,xlab="time in months",ylab="log hazard ratio",
     ylim=c(min(DPresult$loghrl[3,],na.rm=TRUE),
            max(DPresult$loghru[3,],na.rm=TRUE)))
lines(DPresult$predtime,DPresult$loghrl[3,],
      lwd=1,lty=1)
lines(DPresult$predtime,DPresult$loghru[3,],
      lwd=1,lty=1)
dev.off()

pdf("LdataStage4v1HRptest.pdf")
plot(DPresult$predtime,exp(DPresult$loghr.est[3,]),
     main="Hazard ratio for Stage 4 vs Stage 1",
     type="l",lwd=4,xlab="time in months",ylab="hazard ratio",
     ylim=c(min(exp(DPresult$loghr.est[3,]),na.rm=TRUE),
            max(exp(DPresult$loghr.est[3,]),na.rm=TRUE)))
dev.off()

## plot(DPresult$predtime,DPresult$loghr.est[4,],
##      main="log hazard ratio for Age",
##      type="l",lwd=10,xlab="time",ylab="",
##      ylim=c(min(DPresult$loghrl[4,],na.rm=TRUE),
##             max(DPresult$loghru[4,],na.rm=TRUE)))
## lines(DPresult$predtime,DPresult$loghrl[4,],
##       lwd=10,lty=3)
## lines(DPresult$predtime,DPresult$loghru[4,],
##       lwd=10,lty=3)

pdf("LdataStage1and4Survival.pdf")
plot(DPpredict$tpred,DPpredict$Spred[1,],
     main="Estimated Survival for Stage=1 and 4,\n Age=65, Year of Diaganosis=74",
     type="l",lwd=4,xlab="time in months",ylab="survival probability",
     ylim=c(0,1),xlim=c(0,15))
lines(DPpredict$tpred,DPpredict$Spredl[1,],lwd=1,lty=1,type="l")
lines(DPpredict$tpred,DPpredict$Spredu[1,],lwd=1,lty=1,type="l")
lines(DPpredict$tpred,DPpredict$Spred[4,],lwd=4,lty=4,type="l")
lines(DPpredict$tpred,DPpredict$Spredl[4,],lwd=1,lty=4,type="l")
lines(DPpredict$tpred,DPpredict$Spredu[4,],lwd=1,lty=4,type="l")
dev.off()

medianSurv <- function(time,S){
   time[apply(S,1,function(x){findInterval(-0.5,-x)})]
}

pdf("LdataStage1and4MedSurvTime.pdf")
medianTime1 <- medianSurv(DPpredict$tpred,DPpredict$S[,,1])
medianTime4 <- medianSurv(DPpredict$tpred,DPpredict$S[,,4])
postmed1 <- density(medianTime1)
postmed4 <- density(medianTime4)
plot(postmed4,main="Posterior of Median Survival, Stage=1 and 4",
     type="l",lwd=4,xlab="median survival time in months",ylab="density",
     xlim=c(0,15))
lines(postmed1,type="l",lwd=4,lty=4)     
## hist(medianTime1,breaks=40,xlim=c(4,15),probability = TRUE,
##      main="Histogram of median survival time")
abline(v=quantile(medianTime4,0.025),lwd=2,lty=1)
abline(v=quantile(medianTime4,0.975),lwd=2,lty=1)
abline(v=quantile(medianTime1,0.025),lwd=2,lty=4)
abline(v=quantile(medianTime1,0.975),lwd=2,lty=4)
legend(8,0.4,c("Stage 4","Stage 1"),lty=c(1,4),lwd=4)
dev.off()
