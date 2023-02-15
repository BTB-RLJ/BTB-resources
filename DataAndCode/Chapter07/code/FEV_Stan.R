# Loads the R packages with appropriate functions into the R session
library(rstan)
library(coda)
library(ggmcmc)

options(mc.cores = parallel::detectCores()) # Allows Stan to use multiple cores

# Data
FEVdata = read.table("../data/FEVdata.txt", header=TRUE)

mtilde <- c(2.8,3,4,3.3) # Specify prior mean vector
D <- diag(c(0.04,0.04,0.04,0.09)) # and covariance matrix for mtilde
Xtilde <- matrix(c(1,11,0,0,  # Specify X tilde matrix
                   1,13,1,13,
                   1,16,0,0,
                   1,18,1,18),4,4,byrow=T)
Xtildeinv <- solve(Xtilde)   # Invert X tilde
# Get prior mean vector and precision matrix for beta
beta0 <- c(t(Xtildeinv %*% mtilde))
C0 <- Xtildeinv %*% D %*% t(Xtildeinv)
C0inv <- solve(C0)
a <- 1.73   # a and b are the hyperparameters of prior on tau
b <- 0.78
n <- length(FEVdata$FEV)   # Get the sample size, n
p <-  dim(Xtildeinv)[1]  # Number of regression coefficients
# Create a list of all the inputs appearing in the BUGS code
data <- list('n'=n, 'p'=p, 'beta0'=beta0, 'C0'=C0,
             'a'=a, 'b'=b, 'FEV'=FEVdata$FEV, 'Age'=FEVdata$Age, 'Smoke'=FEVdata$Smoke)

# Identify all objects one might want to monitor in MCMC program
parameters <- c('beta','tau','meanFEVs','meanFEVns',
                 'RM','MD','FEV20s','FEV20ns')

# Specify initial values for all stochastic nodes
inits <- function(){
	list(beta=rnorm(4, c(2.5,0,0,0), 1), tau=rgamma(1, 1, 1))
}

stan_fit <- stan(file = 'FEV.stan', data = data, chains = 3, iter=2000, warmup=1000,
                 pars = c('beta', 'sigma', 'FEV20s', 'FEV20ns', 'RM', 'MD'),
                 init = inits)
                 
stan_df <- ggs(stan_fit)

stan_df %>%
    group_by(Parameter) %>%
    summarize(mean = mean(value),
              sd = sd(value),
              `2.5%` = quantile(value, .025),
              median = median(value),
              `97.5%` = quantile(value, .975))
