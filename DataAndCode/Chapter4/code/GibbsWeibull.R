# 9/8/20 Tested; works. 
# 8/29/20
# Notation now consistent with modified introduction of uniform-Pareto
# in Section 4.3.1 and in Chapter 11
require("HI")
logdensity <- function(alpha,phi,y,n,logprody,lambda)
  (n*log(alpha)+alpha*logprody-lambda*sum(exp(alpha*log(y))))
gibbsWeibull <- function(
  y=rexp(500)^0.5,  # default data from Weib(alpha=2,lambda=1)
  a=1.0,            # shape parameter in prior for lambda
  b=1.0,            # rate parameter in prior for lambda
  c=1,              # Pareto lower bound parameter
  d=2,              # Pareto shape parameter
  burnin=100,       # initial iterations to discard
  thin=1,           # thining factor; no thining by default
  M=1000,           # number of samples returned
  lambdainit=NULL,  # initial value for lambda
  alphainit=NULL,   # initial value for Weibull shape parameter
  phiinit=NULL      # initial value for hierarchical parameter
  )
{
  if(length(lambdainit)==0) lambdainit <- a/b
  if(length(phiinit)==0) phiinit <- c*d/(d-1)
  if(length(alphainit)==0) alphainit <- phiinit/2
  lambdas <- rep(0,M);  alphas <- rep(1,M)
  lambda <- lambdainit; alpha <- alphainit; phi <- phiinit
  n <- length(y);  logprody <- sum(log(y))
  a1 <- a+n
  d1 <- d+1
  for(g in (1-burnin):(M*thin)){ # begin Gibbs loop
    b1 <- b+sum(exp(alpha*log(y)))
    lambda <- rgamma(1,shape=a1,rate=b1)
    phi <- max(c,alpha)/runif(1)^(1/d1)
    alpha <- arms(alpha,logdensity,
                  function(x,...) ((x>0)*(x<phi)),1,phi,
                  y,n,logprody,lambda)
    if(g>0 & (g%%thin)==0){      # save posterior sample
      lambdas[g] <- lambda;   alphas[g] <- alpha
    }
  }                              # end Gibbs loop
  return(list(lambda=lambdas,alpha=alphas))
}

samples <- gibbsWeibull()
summary(samples$lambda)
summary(samples$alpha)

## > logdensity <- function(alpha,phi,y,n,logprody,lambda)
## +   (n*log(alpha)+alpha*logprody-lambda*sum(exp(alpha*log(y))))
## > gibbsWeibull <- function(
## +   y=rexp(500)^0.5,  # default data from Weib(alpha=2,lambda=1)
## +   a=1.0,            # shape parameter in prior for lambda
## +   b=1.0,            # rate parameter in prior for lambda
## +   c=1,              # Pareto lower bound parameter
## +   d=2,              # Pareto shap .... [TRUNCATED] 
## > samples <- gibbsWeibull()
## summary(samples$lambda)
## summary(samples$alpha)
## > 
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##  0.8671  0.9562  0.9899  0.9903  1.0220  1.1488 
## >  
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##   1.754   1.971   2.021   2.020   2.070   2.219 
## 
## > plot(density(samples$lambda))
## 
## > plot(density(samples$alpha))
## >
