library(rstan)
library(rjags)
library(coda)
library(ggmcmc)


data <- list(n=c(19,15,50),
     low=c(91,46,95,60,33,410,105,43,189,1097,54,178,114,137,
           233,101,25,70,357),
     normal=c(370,267,99,157,75,1281,48,298,268,62,804,
              430,171,694,404),
     high=c(75,52,1378,555,331,231,472,263,120,46,650,349,251,
            492,759,96,627,171,1584,69,368,509,486,354,351,
            839,88,162,1041,383,234,1130,503,244,606,457,460,
            283,767,576,628,239,583,428,452,723,201,406,422,243))
# Initial values
inits <- list(mu=c(4.87,5.39,6.4),tau=1)


jags_model <- jags.model(file = "DiasorinANOVA.txt", data = data,
                         n.chains = 3, inits = inits)

ndraws <- 1000
burnin <- 1000
update(jags_model, burnin)
jags_output <- coda.samples(jags_model,
                            variable.names = c('med'),
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
stan_fit <- stan(file = 'DiasorinANOVA.stan', data = data, chains = 3)
stan_df <- ggs(stan_fit, family = "med")
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


