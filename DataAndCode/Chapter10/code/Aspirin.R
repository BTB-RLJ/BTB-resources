library(rstan)
library(rjags)
library(coda)
library(ggmcmc)

aspirin  = read.table("aspirin_data.txt", header=F, 
                      col.names=c("ID", "X", "Y"))
PreMinusPost = aspirin$X - aspirin$Y

data <- list(N=nrow(aspirin), T0=aspirin$X, meanT0 = mean(aspirin$X),
             Y=PreMinusPost, nu = 4)
inits <- list(
    b0 = rnorm(1, 0, 1), b1 = rnorm(1, 0, 1),
    Precision = rgamma(1, 1, 1),
    lambda = runif(12, 3.5, 4.5))
jags_model <- jags.model(file = "Aspirin.txt", data = data, inits = inits,
                         n.chains = 3)

ndraws <- 5000
burnin <- 1000
update(jags_model, burnin)
jags_output <- coda.samples(jags_model,
                            variable.names = c('b0', 'b1', 'StDev'),
                            n.iter=ndraws)
jags_df <- ggs(jags_output)
jags_df %>%
    group_by(Parameter) %>%
    summarize(median = median(value),
              mean = mean(value),
              `2.5%` = quantile(value, .025),
              `97.5%` = quantile(value, .975))

### Now stan
inits$sigma2 <- 1 / inits$Precision
stan_fit <- stan(file = 'Aspirin.stan', data = data, chains = 3,
                 init = rep(list(inits), 3), iter = 10000)
stan_df <- ggs(stan_fit) %>%
    filter(Parameter %in% c('b0', 'b1', 'StDev'))
stan_df %>%
    group_by(Parameter) %>%
    summarize(mean = mean(value),
              `2.5%` = quantile(value, .025),
              median = median(value),
              `97.5%` = quantile(value, .975))

df <- bind_rows(mutate(jags_df, sampler = "JAGS"),
                mutate(stan_df, sampler = "STAN"))
ggplot(df, aes(x = value, color = sampler)) +
    geom_line(stat = "density") +
    facet_wrap(~Parameter, scales = "free")



