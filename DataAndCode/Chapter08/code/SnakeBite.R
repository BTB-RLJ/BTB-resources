library(rstan)
library(rjags)
library(coda)
library(ggmcmc)


data <- list(y = c(8,5,3,0), n = c(8,8,8,8), sDose = c(6,2,0,-1))

test <- glm(cbind(data$y, data$n - data$y) ~ data$sDose, family = "binomial")


jags_model <- jags.model(file = "SnakeBite.jags", data = data,
                         n.chains = 3)

ndraws <- 5000
burnin <- 1000
update(jags_model, burnin)
jags_output <- coda.samples(jags_model,
                            variable.names = c('p', 'ED50'),
                            n.iter=ndraws)
jags_df <- ggs(jags_output)
jags_df %>%
    group_by(Parameter) %>%
    summarize(`2.5%` = quantile(value, .025),
              median = median(value),
              `97.5%` = quantile(value, .975))

### Now stan
stan_fit <- stan(file = 'SnakeBite.stan', data = data, chains = 3)
stan_df <- ggs(stan_fit, family = "p")
stan_df %>%
    group_by(Parameter) %>%
    summarize(`2.5%` = quantile(value, .025),
              median = median(value),
              `97.5%` = quantile(value, .975))

df <- bind_rows(mutate(jags_df, sampler = "JAGS"),
                mutate(stan_df, sampler = "STAN"))
ggplot(df, aes(x = value, color = sampler)) +
    geom_line(stat = "density") +
    facet_wrap(~Parameter, scales = "free")


