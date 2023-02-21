library(rstan)
library(rjags)
library(coda)
library(ggmcmc)

df <- read.table("MarijuanaData.txt")
colnames(df) <- c("age", "n", "l", "u")

data <- list(n = df$n, r = df$n, l = df$l, u = df$u)
inits = list(mu=0, tau=1)

jags_model <- jags.model(file = "Marijuana.txt", data = data,
                         n.chains = 3, inits = inits)

ndraws <- 1000
burnin <- 1000
update(jags_model, burnin)
jags_output <- coda.samples(jags_model,
                            variable.names = c('mu', 'sigma'),
                            n.iter=ndraws)
jags_df <- ggs(jags_output)
jags_df %>%
    group_by(Parameter) %>%
    summarize(mean = mean(value),
              median = median(value),
              `2.5%` = quantile(value, .025),
              `97.5%` = quantile(value, .975))

### Now stan
stan_fit <- stan(file = 'Marijuana.stan', data = data, chains = 3,
                 init = rep(list(inits), 3))
stan_df <- ggs(stan_fit) %>%
    filter(Parameter %in% c('mu', 'sigma'))
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



