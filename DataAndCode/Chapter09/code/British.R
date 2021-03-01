library(rstan)
library(rjags)
library(coda)
library(ggmcmc)

data <- list(A=c(1,2,3,4,5,1,2,3,4,5),
             S=c(1,1,1,1,1,0,0,0,0,0),
             y=c(32,104,206,186,102, 2,12,28,28,31),
             M=c(5.2407,4.3248,2.8612,1.2663,.5317,1.8790,1.0673,.5710,.2585,.1462) )

jags_inits <- list(lam = c( 50, 80))

jags_model <- jags.model(file = "British.txt", data = data,
                         n.chains = 3, inits = jags_inits)

ndraws <- 1000
burnin <- 1000
update(jags_model, burnin)
jags_output <- coda.samples(jags_model,
                            variable.names = c('RR'),
                            n.iter=ndraws)
jags_df <- ggs(jags_output)
jags_df %>%
    group_by(Parameter) %>%
    summarize(mean = mean(value),
              median = median(value),
              `2.5%` = quantile(value, .025),
              `97.5%` = quantile(value, .975))

### Now stan
b10  <- log(50)
sig10 <- (log(100) - b10)/1.645
tau10 <- 1/ sig10 ^ 2
b20 <-  log(80)
sig20 <- (log(150) - b20)/1.645
tau20 <- 1/ sig20 ^ 2
data <- c(data, b10 = b10, tau10 = tau10, b20 = b20, tau20 = tau20)
stan_fit <- stan(file = 'British.stan', data = data,
                 chains = 3, init = list(jags_inits, jags_inits, jags_inits))
stan_df <-
    ggs(stan_fit, family = 'RR') 
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






