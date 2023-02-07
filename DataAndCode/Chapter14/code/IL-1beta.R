library(rstan)
library(rjags)
library(coda)
library(ggmcmc)

options(mc.cores = parallel::detectCores()) # Allows Stan to use multiple cores

df <- read.table('../data/IL-1betaData.txt', header=TRUE)
# Data
data <- list(logy = log(df$y),
             N = nrow(df),
             ID = df$ID,
             nID = length(unique(df$ID)))
             
inits = function() {
	list(beta1=rnorm(1, 0, 1), rhow=runif(1, 0.1, 0.9), sig=runif(1, 0.5, 1.5), sigg=runif(1, 0.5, 1.5), sigw=runif(1, 0.5, 1.5))
}

jags_model <- jags.model(file = "IL-1beta.jags", data = data,
                         n.chains = 3, inits = inits)

ndraws <- 8000
burnin <- 5000
update(jags_model, burnin)

jags_output <- coda.samples(jags_model,
                            variable.names = c('yf', 'med', 'ten[2]', 'ninety[2]', 
                            		'sig', 'sigg', 'rho', 'beta1'),
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
stan_fit <- stan(file = 'IL-1beta.stan', data = data, chains = 3, iter=8000, warmup=5000,
                                  init = inits)

stan_df <- ggs(stan_fit) %>%
	filter(Parameter %in% c('yf', 'med', 'ten[2]', 'ninety[2]', 'sig', 'sigg', 'rho', 'beta1'))
#    filter(Parameter %in% c('med', 'sig', 'rho'))

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


