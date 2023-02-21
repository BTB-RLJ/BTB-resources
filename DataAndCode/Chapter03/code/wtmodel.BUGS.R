library(R2WinBUGS)
# library(R2OpenBUGS) # If using OpenBUGS
library(coda)
library(ggmcmc)

VetBPdata = read.csv("../data/VetBP/VetBP.csv", sep="/"), header=T)

ok = !is.na(VetBPdata$Weight_00)
y = VetBPdata$Weight_00[ok]
n = length(y)
alpha0 <- 0.1; lambda0 <- 0.1;
mu0 <- 150; tau0 <- 0.001
data = list("n", "y", "mu0", "tau0", "alpha0", "lambda0")

inits = function() {
  list(mu = rnorm(1, mean=120, sd=5), tau = runif(1, min=0.1, max=10))
}

parameters = c("mu", "sigma")

# Run with WinBUGS
VetBP.sim = bugs(data, inits, parameters, "wtmodel.txt", 
                n.chains=3, n.iter=10000, debug=F, useWINE=F)
