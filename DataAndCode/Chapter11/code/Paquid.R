library(prodlim)
library(riskRegression)
library(DPWeibull)
data(Paquid)
setwd("C:/Users/ys8wp/Desktop/")
CRresult<-dpweib(Hist(time, status)~1,data=Paquid,
                  predtime = seq(from=0,to=max(Paquid$time),length=200))

save(CRresult,file="CRresult.RData")

pdf("competingRisks.pdf",height=20,width=30,pointsize=30)
plot(CRresult$predtime,CRresult$CIF1.est,
     main="Cumulative Incidence Function",
     type="l",lwd=10,xlab="time",ylab="",col="red",
     ylim=c(0,max(c(CRresult$CIF1u,CRresult$CIF2u),na.rm=TRUE)))
lines(CRresult$predtime,CRresult$CIF1u,
      lwd=10,lty=2,col="red")
lines(CRresult$predtime,CRresult$CIF1l,
      lwd=10,lty=2,col="red")
lines(CRresult$predtime,CRresult$CIF2u,
      lwd=10,lty=2,col="blue")
lines(CRresult$predtime,CRresult$CIF2l,
      lwd=10,lty=2,col="blue")
lines(CRresult$predtime,CRresult$CIF2.est,
      lwd=10,col="blue")
legend("topleft",c("Cause 1","Cause 2"),col=c("red","blue"),lty=1,lwd=10)

plot(CRresult$predtime,CRresult$d1.est,
     main="Subdistribution Density Function",
     type="l",lwd=10,xlab="time",ylab="",col="red",
     ylim=c(0,max(c(CRresult$d1u,CRresult$d2u),na.rm=TRUE)))
lines(CRresult$predtime,CRresult$d1u,
      lwd=10,lty=2,col="red")
lines(CRresult$predtime,CRresult$d1l,
      lwd=10,lty=2,col="red")
lines(CRresult$predtime,CRresult$d2u,
      lwd=10,lty=2,col="blue")
lines(CRresult$predtime,CRresult$d2l,
      lwd=10,lty=2,col="blue")
lines(CRresult$predtime,CRresult$d2.est,
      lwd=10,col="blue")

plot(CRresult$predtime,CRresult$h1.est,
     main="Subdistribution Hazard Function",
     type="l",lwd=10,xlab="time",ylab="",col="red",
     ylim=c(0,max(c(CRresult$h1u,CRresult$h2u),na.rm=TRUE)))
lines(CRresult$predtime,CRresult$h1u,
      lwd=10,lty=2,col="red")
lines(CRresult$predtime,CRresult$h1l,
      lwd=10,lty=2,col="red")
lines(CRresult$predtime,CRresult$h2u,
      lwd=10,lty=2,col="blue")
lines(CRresult$predtime,CRresult$h2l,
      lwd=10,lty=2,col="blue")
lines(CRresult$predtime,CRresult$h2.est,
      lwd=10,col="blue")
dev.off()