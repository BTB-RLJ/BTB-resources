library(rstan)
library(rjags)
library(coda)
library(ggmcmc)

options(mc.cores = parallel::detectCores()) # Allows Stan to use multiple cores

### Read data
VetBPdata = read.csv('../data/VetBP.csv', header=T)

x <- VetBPdata$Age - 60
y <- VetBPdata$DBP_00
n <- length(x)
mubeta0 <- 75
taubeta0 <- 0.01
mubeta1 <- 0
taubeta1 <- 4
a <- 5
b <- 400
data <- list('n'=n, 'y'=y, 'x'=x, 'mubeta0'=mubeta0, 'taubeta0'=taubeta0,
                'mubeta1'=mubeta1, 'taubeta1'=taubeta1, 'a'=a, 'b'=b)
                
inits <- function(){
   list(beta0=runif(1,50,100), beta1=runif(1,-1,1), tau=runif(1,0.001,0.03))
}

parameters <- c('beta0','beta1','tau', 'sigma')

jags_model <- jags.model(file = 'dbpagemodel.jags', data = data,
                         n.chains = 3, inits = inits)

ndraws <- 2000
burnin <- 2000

update(jags_model, burnin)

jags_output <- coda.samples(jags_model,
                            variable.names = c('beta0','beta1','tau', 'sigma'),
                            n.iter=ndraws)
                            
jags_df <- ggs(jags_output)

### Check convergence
ggs_traceplot(jags_df) + facet_wrap(~Parameter, scales = 'free_y')

jags_df %>%
    group_by(Parameter) %>%
    summarize(mean = mean(value),
              sd = sd(value),
              `2.5%` = quantile(value, .025),
              median = median(value),
              `97.5%` = quantile(value, .975))

### Now stan
stan_fit <- stan(file = 'dbpagemodel.stan', data = data, chains = 3,
					pars=c('beta0','beta1','tau', 'sigma'),
					init=inits, iter=4000, warmup=2000)

stan_df <- ggs(stan_fit)

### Check convergence
ggs_traceplot(stan_df) + facet_wrap(~Parameter, scales = 'free_y')

stan_df %>%
    group_by(Parameter) %>%
    summarize(mean = mean(value),
              sd = sd(value),
              `2.5%` = quantile(value, .025),
              median = median(value),
              `97.5%` = quantile(value, .975))

df <- bind_rows(mutate(jags_df, sampler = 'JAGS'),
                mutate(stan_df, sampler = 'STAN'))
ggplot(df, aes(x = value, color = sampler)) +
    geom_line(stat = 'density') +
    facet_wrap(~Parameter, scales = 'free')


