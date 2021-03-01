library(rstan)
library(rjags)
library(coda)
library(ggmcmc)


# Data
trauma <- read.table("http://people.oregonstate.edu/~calverta/BIDA/Chapter08/trauma300.txt",
                     header = TRUE)

Xtilde <- c(1,25,7.84,60,0,0,
            1,25,3.34,10,0,0,
            1,41,3.34,60,1,60,
            1,41,7.84,10,1,10,
            1,33,5.74,35,0,0,
            1,33,5.74,35,1,35)
Xtilde <- matrix(Xtilde, nrow=6, ncol=6, byrow=TRUE)
Xtinv <- solve(Xtilde)

data <- list(death = trauma$death, ISS = trauma$ISS,
             TI = trauma$TI, RTS = trauma$RTS, AGE = trauma$AGE,
             n=300, a=c(1.1,3,5.9,1.3,1.1,1.5), b=c(8.5,11,1.7,12,4.9,5.5),
             Xtinv = Xtinv)

jags_model <- jags.model(file = "Trauma.txt", data = data,
                         n.chains = 3)

ndraws <- 1000
burnin <- 1000
update(jags_model, burnin)
jags_output <- coda.samples(jags_model,
                            variable.names = c('beta'),
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
stan_fit <- stan(file = 'Trauma.stan', data = data, chains = 3)
stan_df <- ggs(stan_fit, family = "beta")

df <- bind_rows(mutate(jags_df, sampler = "JAGS"),
                mutate(stan_df, sampler = "STAN"))
ggplot(df, aes(x = value, color = sampler)) +
    geom_line(stat = "density") +
    facet_wrap(~Parameter, scales = "free")


