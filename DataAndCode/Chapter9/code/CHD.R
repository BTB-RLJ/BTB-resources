library(rstan)
library(rjags)
library(coda)
library(ggmcmc)

data <- list(M = c(42.259, 61.451), y = c(823, 650), Sex = c(0,1))

jags_inits <- list(b = c(0,0)) 

jags_model <- jags.model(file = "CHD.txt", data = data,
                         n.chains = 3, inits = jags_inits)

ndraws <- 1000
burnin <- 1000
update(jags_model, burnin)
jags_output <- coda.samples(jags_model,
                            variable.names = c('lam', 'RR'),
                            n.iter=ndraws)
jags_df <- ggs(jags_output)
jags_df %>%
    group_by(Parameter) %>%
    summarize(mean = mean(value),
              median = median(value),
              `2.5%` = quantile(value, .025),
              `97.5%` = quantile(value, .975))

### Now stan
stan_fit <- stan(file = 'CHD.stan', data = data,
                 chains = 3, init = 0)
stan_df <-
    ggs(stan_fit) %>%
    filter(Parameter %in% c('lam[1]', 'lam[2]', 'RR'))
stan_df %>%
    group_by(Parameter) %>%
    summarize(mean = mean(value),
              median = median(value),
              `2.5%` = quantile(value, .025),
              `97.5%` = quantile(value, .975))

df <- bind_rows(mutate(jags_df, sampler = "JAGS"),
                mutate(stan_df, sampler = "STAN"))
ggplot(df, aes(x = value, color = sampler)) +
    geom_line(stat = "density") +
    facet_wrap(~Parameter, scales = "free")






