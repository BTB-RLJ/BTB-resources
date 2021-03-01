library(rstan)
library(rjags)
library(coda)
library(ggmcmc)


# Data
data <- list(High=c(1,0), y = c(126 ,53 ),
             n=c(477 ,734 ), a=c(2.044,2.044), b=c(13.006,13.006),
             Xtinv = matrix(c(0, 1,
                              1, -1),
                            nrow  = 2, byrow = TRUE))
jags_inits <- list(tildetheta=c(0.5,0.5))

jags_model <- jags.model(file = "LymphedemaBF.txt", data = data,
                         n.chains = 3, inits = jags_inits)

ndraws <- 1000
burnin <- 1000
update(jags_model, burnin)
jags_output <- coda.samples(jags_model,
                            variable.names = c('w'),
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
stan_fit <- stan(file = 'LymphedemaBF.stan', data = data, chains = 3,
                 init = rep(list(jags_inits), 3))
stan_df <- ggs(stan_fit, family = "w")
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


