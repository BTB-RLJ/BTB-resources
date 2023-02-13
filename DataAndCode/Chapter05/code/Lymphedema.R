library(rstan)
library(rjags)
library(coda)
library(ggmcmc)

options(mc.cores = parallel::detectCores()) # Allows Stan to use multiple cores

### Data for lymphedema example
data = list(n=c(477, 734), y=c(126, 53), a=c(2.044, 13.006), b=c(2.044, 13.006))

# Initial values
inits <- function(){
	list(theta=rbeta(2, 4, 4))
}

jags_model <- jags.model(file = "Lymphedema.jags", data = data,
                         n.chains = 3, inits = inits)

ndraws <- 2000
burnin <- 5000

update(jags_model, burnin)

jags_output <- coda.samples(jags_model,
                            variable.names = c('theta', 'RR', 'OR', 'RD', 'test'),
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
'
### Now stan
stan_fit <- stan(file = 'Lymphedema.stan', data = data, chains = 3,
					pars=c('theta', 'RR', 'OR', 'RD', 'test'),
					init=inits, iter=7000, warmup=5000)

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


