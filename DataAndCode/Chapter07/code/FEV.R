 library(R2WinBUGS) # Loads the R package with appropriate
                    # functions into the R session
 tilde m <- c(2.8,3,4,3.3) # Specify prior mean vector
 D <- diag(c(0.04,0.04,0.04,0.09)) # and cov matrix for mtilde
 Xtilde <- matrix(c(1,11,0,0,  # Specify tilde X matrix
                    1,13,1,13,
                    1,16,0,0,
                    1,18,1,18),4,4,byrow=T)
 Xtildeinv <- solve(Xtilde)   # Invert tilde X
 # Get prior mean vector and precision matrix for beta
 beta0 <- c(t(Xtildeinv %*% tilde m))
 C0 <- Xtildeinv %*% D %*% t(Xtildeinv)
 C0inv <- solve(C0)
 a <- 1.73   # a and b are the hyperparameters of prior on tau
 b <- 0.78
 n <- length(FEV)   # Get the sample size, n
 p <-  dim(Xtildeinv)[1]  # Number of regression coefficients
 # Create a list of all the inputs appearing in the BUGS code
 FEVdataBUGS <- list("n","r","beta0","C0inv",
                     "a","b","FEV","Age","Smoke")
 # Identify all objects to be monitored in WinBUGS
 parameters <- c("beta","tau","meanFEVs","meanFEVns",
                 "RM","MD","FEV20s","FEV20ns")
 # Specify initial values for all stochastic nodes
 inits <- list(list(tau=1,beta=c(0,0,0,0),FEV20s=2, FEV20ns=2))
 FEV.fit <- bugs(FEVdataBUGS, inits, parameters,"FEV.bug",
            n.chains=1, n.iter=60000, n.thin=1, n.burnin=10000)
 print(FEV.fit,digits=3)
 attach.bugs(FEV.fit)

