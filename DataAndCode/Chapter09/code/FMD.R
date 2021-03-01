library(rstan)
library(rjags)
library(coda)
library(ggmcmc)
df <- read.table("http://people.oregonstate.edu/~calverta/BIDA/Chapter11/FMDData.txt",
                   header = TRUE)

data <- list(Size = df$Cattle,
             East = df$EasternTurkey,
             FMD = df$FMD1998)


jags_model <- jags.model(file = "FMD.txt", data = data,
                         n.chains = 3)

ndraws <- 1000
burnin <- 1000
update(jags_model, burnin)
jags_output <- coda.samples(jags_model,
                            variable.names = c('b'),
                            n.iter=ndraws)
jags_df <- ggs(jags_output)
jags_df %>%
    group_by(Parameter) %>%
    summarize(mean = mean(value),
              median = median(value),
              `2.5%` = quantile(value, .025),
              `97.5%` = quantile(value, .975))

### Now stan
stan_fit <- stan(file = 'FMD.stan', data = data,
                 chains = 3)
stan_df <- ggs(stan_fit, family = '^b')
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






