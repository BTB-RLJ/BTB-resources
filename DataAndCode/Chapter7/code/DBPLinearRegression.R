library(rstan)
library(rjags)
library(coda)
library(ggmcmc)

vetBP <- read.table("VetBP.txt")

data <- list(DBP = vetBP$V5,
             Bdp = vetBP$V4,
             Wt = vetBP$V3,
             Ht = vetBP$V2,
             Trt = vetBP$V1,
             N = nrow(vetBP))
jags_inits <- list(b = c(0,0,0,0,0), gamma =0) 

jags_model <- jags.model(file = "DBPLinearRegression.jag", data = data,
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
              sd = sd(value),
              `2.5%` = quantile(value, .025),
              `97.5%` = quantile(value, .975))

### Now stan
stan_fit <- stan(file = 'DBPLinearRegression.stan', data = data,
                 chains = 3)
stan_df <- ggs(stan_fit, family = "b")
stan_df %>%
    group_by(Parameter) %>%
    summarize(mean = mean(value),
              sd = sd(value),
              `2.5%` = quantile(value, .025),
              `97.5%` = quantile(value, .975))

df <- bind_rows(mutate(jags_df, sampler = "JAGS"),
                mutate(stan_df, sampler = "STAN"))
ggplot(df, aes(x = value, color = sampler)) +
    geom_line(stat = "density") +
    facet_wrap(~Parameter, scales = "free")


