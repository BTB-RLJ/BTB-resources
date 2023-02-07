library(rstan)
library(rjags)
library(coda)
library(ggmcmc)

options(mc.cores = parallel::detectCores())

df <- read.table('../data/CowAbortionData.txt', sep = "\t", header=T, stringsAsFactors = F, nrows = 13145)
colnames(df) <- c('herd', 'DO', 'GR', 'y')


# Data
data <- list(herd = df$herd, 
             DO = df$DO,
             GR = df$GR,
             y = df$y,
             N = nrow(df))

inits = function(){
	list(beta=rnorm(4, 0, 1), sigmag=runif(1, 0.5, 2.0))
}

jags_model <- jags.model(file = "CowAbortion.jags", data = data,
                         n.chains = 3, inits = inits)

ndraws <- 3000
burnin <- 1000
update(jags_model, burnin)
jags_output <- coda.samples(jags_model,
                            variable.names = c('beta', 'sigmag'),
                            n.iter=ndraws)
jags_df <- ggs(jags_output)

### Show convergence
ggs_traceplot(jags_df) + facet_wrap(~Parameter, scales = "free_y")

jags_df %>%
    group_by(Parameter) %>%
    summarize(mean = mean(value),
              sd = sd(value),
              `2.5%` = quantile(value, .025),
              median = median(value),
              `97.5%` = quantile(value, .975))

### Now stan
stan_fit <- stan(file = 'CowAbortion.stan', data = data, chains = 3,
                 init = inits,
                 pars=c('gamma', 'theta', 'prob'), include=FALSE) # Do not save "gamma" or "p" vectors
stan_df <- ggs(stan_fit) 

### Show convergence
ggs_traceplot(stan_df) + facet_wrap(~Parameter, scales = "free_y")

stan_df %>%
    group_by(Parameter) %>%
    summarize(mean = mean(value),
              sd = sd(value),
              `2.5%` = quantile(value, .025),
              median = median(value),
              `97.5%` = quantile(value, .975))

### Compare results for JAGS and Stan
df <- bind_rows(mutate(jags_df, sampler = "JAGS"),
                mutate(stan_df, sampler = "STAN"))
ggplot(df, aes(x = value, color = sampler)) +
    geom_line(stat = "density") +
    facet_wrap(~Parameter, scales = "free")


