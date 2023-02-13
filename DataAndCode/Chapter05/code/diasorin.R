library(rstan)
library(rjags)
library(coda)
library(ggmcmc)

options(mc.cores = parallel::detectCores()) # Allows Stan to use multiple cores

data <- list(nL=19, nN=15, c=1,d=1,e=1,f=1, bL = 1, bN =1,
 LowTurnover = c(91, 46,  95,  60, 33, 410, 105, 43,
 189,1097, 54,178, 114, 137, 233, 101, 25,70,357),
 NormalTurnover = c(370, 267, 99,157, 75,1281, 48, 298, 268,
 62,804,430,171,694,404))

# Initial values
inits <- function(){
	list(muL=rnorm(1, 4.87, 3), muN=rnorm(1, 5.39, 1), tauL=rgamma(1, 1, 1), tauN=rgamma(1, 1, 1))
}

jags_model <- jags.model(file = "diasorin.jags", data = data,
                         n.chains = 3, inits = inits)

ndraws <- 2000
burnin <- 2000

update(jags_model, burnin)

jags_output <- coda.samples(jags_model,
                            variable.names = c('medL', 'medN', 'test', 'AUC'),
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
stan_fit <- stan(file = 'diasorin.stan', data = data, chains = 3,
					pars=c('medL', 'medN', 'test', 'AUC'),
					init=inits, iter=4000, warmup=2000)

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


