library(rstan)
library(rjags)
library(coda)
library(ggmcmc)

options(mc.cores = parallel::detectCores()) # Allows Stan to use multiple cores

y = c(815, 115, 208, 327)
n = sum(y)
aSe = bSe = 1
aSp = bSp = 1
aPi = bPi = 1

data <- list(n=n, y = y, aSe=aSe, bSe=bSe, aSp=aSp, bSp=bSp, aPi=aPi, bPi=bPi)

inits <- function(){
	list(Se = runif(1, 0.25, 0.75), Sp = runif(1, 0.25, 0.75),  pi = runif(1, 0.25, 0.75))
}

jags_model <- jags.model(file = "CAD.jags", data = data, inits = inits,
                         n.chains = 3)

ndraws <- 5000
burnin <- 1000

update(jags_model, burnin)

jags_output <- coda.samples(jags_model,
                            variable.names = c('pi', 'Se', 'Sp'),
                            n.iter=ndraws)
                            
jags_df <- ggs(jags_output)

### Can check convergence
ggs_traceplot(jags_df) + facet_wrap(~Parameter, scales = "free_y")

jags_df %>%
    group_by(Parameter) %>%
    summarize(mean = mean(value),
              sd = sd(value),
              `2.5%` = quantile(value, .025),
              median = median(value),
              `97.5%` = quantile(value, .975))

### Now stan
stan_fit <- stan(file = 'CAD.stan', data = data, chains = 3,
                 init = inits)
                 
stan_df <- ggs(stan_fit) %>%
    filter(Parameter %in% c('pi', 'Se', 'Sp'))

### Can check convergence
ggs_traceplot(stan_df) + facet_wrap(~Parameter, scales = "free_y")

stan_df %>%
    group_by(Parameter) %>%
    summarize(mean = mean(value),
              sd = sd(value),
              `2.5%` = quantile(value, .025),
              median = median(value),
              `97.5%` = quantile(value, .975))

df <- bind_rows(mutate(jags_df, sampler = "JAGS"),
                mutate(stan_df, sampler = "STAN"))
ggplot(df, aes(x = value, color = sampler)) +
    geom_line(stat = "density") +
    facet_wrap(~Parameter, scales = "free")



