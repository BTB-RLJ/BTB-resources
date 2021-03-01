# Prepare data; read from file named VetBP.csv
ds <- read.csv(VetBP.csv,header=TRUE)
attach(ds)
x <- Age-60
y <- DBP_00
n <- length(x)
mubeta1 <- 75; taubeta1 <- 0.01
mubeta2 <- 0; taubeta2 <- 4
a <- 5; b <- 400
data <- 
list("n","y","x","mubeta0","taubeta0",
                "mubeta1", "taubeta1", "a","b")
inits <- function(){
list(beta0=runif(1,50,100),beta1=runif(1,-1,1), 
                       tau=runif(1,0.001,0.03))
}
parameters <- c("beta0","beta1","tau")

# Run model with OpenBUGS
library("R2WinBUGS")
# BUGS model is in file dbpagemodel.txt
dbponage <- 
     bugs(data,inits,parameters,"dbpagemodel.bug",
     n.chains=3,n.iter=10000,program="openbugs",
     debug=TRUE)

# Convert tau samples to sigma and create an mcmc object
sigma <- 1/sqrt(dbponage$sims.matrix[,3])
mcout <- as.mcmc(cbind(dbponage$sims.matrix[,1:2],sigma))

# Obtain posterior summaries and plots
summary(mcout)
pdf("plots.pdf")
plot(mcout)
autocorr.plot(mcout)
dev.off()

