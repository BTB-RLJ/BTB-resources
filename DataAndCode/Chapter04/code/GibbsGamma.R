# 9/8/20 Tested; works.
# 8/29/20
# Notation now consistent with modified introduction of uniform-Pareto
# in Section 4.3.1 and in Chapter 11

require("HI")
logdensity <-
  function(alpha,phi,n,sumy,logprody,lambda)
  (n*alpha*log(lambda)-n*lgamma(alpha)+alpha*logprody
   -lambda*sumy)

gibbsgamma <- function(
  y=rgamma(500,2,1), # default data from Ga(alpha=2,lambda=1)
  a=1.0,             # shape parameter in prior for lambda
  b=1.0,             # rate parameter in prior for lambda
  c=1,               # Pareto lower bound parameter
  d=2,               # Pareto shape parameter
  burnin=100,        # initial iterations to discard
  thin=1,            # thining factor; no thining by default
  M=1000,            # number of samples returned
  lambdainit=NULL,   # initial value for lambda
  alphainit=NULL,    # initial value for Gamma shape parameter
  phiinit=NULL       # initial value for hierarchical parameter
  )
{
  if(length(lambdainit)==0) lambdainit <- a/b
  if(length(phiinit)==0) phiinit <- c*d/(d-1)
  if(length(alphainit)==0) alphainit <- phiinit/2
  lambdas <- rep(0,M);  alphas <- rep(1,M)
  lambda <- lambdainit; alpha <- alphainit; phi <- phiinit
  n <- length(y);  sumy <- sum(y);  logprody <- sum(log(y))
  d1 <- d+1
  b1 <- b + sumy
  for(g in (1-burnin):(M*thin)){ # begin Gibbs loop
    a1 <- a+n*alpha
    lambda <- rgamma(1,shape=a1,rate=b1)
    phi <- max(c,alpha)/runif(1)^(1/d1)
    alpha <- arms(alpha,logdensity,
                  function(x,...) ((x>0)*(x<phi)),1,phi,
                  n,sumy,logprody,lambda)
    if(g>0 & (g%%thin)==0){      # save posterior sample
      lambdas[g] <- lambda;   alphas[g] <- alpha
    }
  }                              # end Gibbs loop
  return(list(lambda=lambdas,alpha=alphas))
}

samples <- gibbsgamma()
summary(samples$lambda)
summary(samples$alpha)

## > require("HI")
## Loading required package: HI
## > logdensity <-
## +   function(alpha,phi,n,sumy,logprody,lambda)
## +   (n*alpha*log(lambda)-n*lgamma(alpha)+alpha*logprody
## +    -lambda*sumy)
## > 
## + . + 
## > gibbsgamma <- function(
## +   y=rgamma(500,2,1), # default data from Ga(alpha=2,lambda=1)
## +   a=1.0,             # shape parameter in prior for lambda
## +   b=1.0,             # rate parameter in prior for lambda
## +   c=1,               # Pareto lower bound parameter
## +   d=2,               # Pareto sha .... [TRUNCATED] 
## > samples <- gibbsgamma()
## summary(samples$lambda)
## summary(samples$alpha)
## > 
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##  0.8431  1.0392  1.0851  1.0874  1.1382  1.3331 
## >  
##   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##   1.682   1.984   2.056   2.064   2.147   2.516 
## > 

