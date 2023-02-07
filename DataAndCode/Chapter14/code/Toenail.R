library(rstan)
library(rjags)
library(coda)
library(ggmcmc)

options(mc.cores = parallel::detectCores())

df <- read.table('../Data/ToenailData.txt', header = TRUE)

# Data
data <- list(y = df$y,
             trt = df$Trt,
             time = df$Time,
             N = nrow(df),
             MyID = df$ID,
             nID = max(unique(df$ID)))

inits = function(){
	list(beta=rnorm(4, 0, 1), gamma=rnorm(data$nID, 0, 0.5), gamma_sd=runif(1, 0.5, 2.0))
}

jags_model <- jags.model(file = "Toenail.jags", data = data,
                         n.chains = 3, inits = inits)

ndraws <- 4000
burnin <- 6000
update(jags_model, burnin)
jags_output <- coda.samples(jags_model,
                            variable.names = c('beta', 'gamma_sd'),
                            n.iter=ndraws)
jags_df <- ggs(jags_output)
jags_df %>%
    group_by(Parameter) %>%
    summarize(mean = mean(value),
              sd = sd(value),
              `2.5%` = quantile(value, .025),
              median = median(value),
              `97.5%` = quantile(value, .975))

### Now stan

stan_fit <- stan(file = 'Toenail.stan', data = data, chains = 3, iter = 10000, warmup = 6000,
                 init = inits, 
                 pars=c("gamma", "p"), include=FALSE) # Do not save "gamma" or "p" vectors
stan_df = ggs(stan_fit)

#	Can also just extract the beta vector and the gamma standard deviation separately
# stan_beta <- ggs(stan_fit, family='beta')
# stan_gamma_sd <- ggs(stan_fit, family='gamma_sd')
# stan_df = rbind(stan_beta, stan_gamma_sd)

### Show convergence
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
