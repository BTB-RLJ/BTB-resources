library(rstan)
library(rjags)
library(coda)
library(ggmcmc)

options(mc.cores = parallel::detectCores()) # Allows Stan to use multiple cores

### Data for Reye's Syndrome example
n = c(7, 16)
x = c(7, 8)

data = list(n=n, x=x, a=1, b=1, mu=1.1, prec=1/(1.03*1.03))

# Initial values
inits <- function(){
	list(delta=rnorm(1, 1, 0.5))
}

jags_model <- jags.model(file = "Reyes.jags", data = data,
                         n.chains = 3, inits = inits)

ndraws <- 2000
burnin <- 4000

update(jags_model, burnin)

jags_output <- coda.samples(jags_model,
                            variable.names = c('delta', 'ttheta', 'OR'),
                            n.iter=ndraws)
jags_df <- ggs(jags_output)

### Check convergence
ggs_traceplot(jags_df) + facet_wrap(~Parameter, scales = "free_y")

jags_df %>%
    group_by(Parameter) %>%
    summarize(mean = mean(value),
              sd = sd(value),
              `2.5%` = quantile(value, .025),
              median = median(value),
              `97.5%` = quantile(value, .975))

### Now stan
stan_fit <- stan(file = 'Reyes.stan', data = data, chains = 3,
					pars=c('delta', 'ttheta', 'OR'),
					init=inits, iter=6000, warmup=4000)

stan_df <- ggs(stan_fit)

### Check convergence
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


