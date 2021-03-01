library(rstan)
library(rjags)
library(coda)
library(ggmcmc)

data <- list(N = 27,
             x=c(1, 1.5, 1.5, 1.5, 2.5, 4, 5, 5, 7, 8, 8.5, 9, 9.5, 9.5, 10, 
                 12, 12, 13, 13, 14.5, 15.5, 15.5, 16.5, 17, 22.5, 29, 31.5),
             Y = c(1.8, 1.85, 1.87, 1.77, 2.02, 2.27, 2.15, 2.26, 2.47, 2.19, 
                   2.26, 2.4, 2.39, 2.41, 2.5, 2.32, 2.32, 2.43, 2.47, 2.56, 2.65, 
                   2.47, 2.64, 2.56, 2.7, 2.72, 2.57))
             
inits <- list(alpha = 1,
              beta = 1,
              tau = 1,
              gamma = 0.9)

jags_model <- jags.model(file = "Dugongs.txt", data = data,
                         n.chains = 3, inits = inits)

ndraws <- 1000
burnin <- 1000
update(jags_model, burnin)
jags_output <- coda.samples(jags_model,
                            variable.names = c('alpha', 'beta', 'tau', 'gamma'),
                            n.iter=ndraws)
jags_df <- ggs(jags_output)
jags_df %>%
    group_by(Parameter) %>%
    summarize(mean = mean(value),
              median = median(value),
              `2.5%` = quantile(value, .025),
              `97.5%` = quantile(value, .975))

### Now stan
stan_fit <- stan(file = 'Dugongs.stan', data = data, chains = 3,
                 init = rep(list(inits), 3))
stan_df <-
    ggs(stan_fit) %>%
    filter(Parameter %in% c('alpha', 'beta', 'tau', 'gamma'))
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



