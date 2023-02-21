library(rstan)
library(rjags)
library(coda)
library(ggmcmc)


# Data
data <- list(n=c(17,16),
             a1=.001, a2=.001, b1=.001, b2=.001, 
             t1=c(65,156,100,134,16,108,121,4,39,143,56,26,22,1,1,5,65),
             t2=c(56,65,17,7,16,22,3,4,2,3,8,4,3,30,4,43))
inits = list(theta1=0.05,theta2=0.02)

jags_model <- jags.model(file = "LeukemiaModeling.txt", data = data,
                         n.chains = 3, inits = inits)

ndraws <- 1000
burnin <- 1000
update(jags_model, burnin)
jags_output <- coda.samples(jags_model,
                            variable.names = c('theta1', 'theta2'),
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
stan_fit <- stan(file = 'LeukemiaModeling.stan', data = data, chains = 3,
                 init = rep(list(inits), 3))
stan_df <- ggs(stan_fit) %>% 
    filter(Parameter %in% c('theta1', 'theta2'))
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


