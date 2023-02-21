library(rstan)
library(rjags)
library(coda)
library(ggmcmc)

options(mc.cores = parallel::detectCores()) # Allows Stan to use multiple cores

ICU <- read.table("../data/Ch14-ICU-Metaanalysis.txt", header=TRUE)
data <- list(m=14, LOR=ICU$LOR, se=ICU$se)

inits <- list(mu=runif(1,0,1),sigmag=runif(1, 0.5, 1.5))

jags_model <- jags.model(file = "Ch14-ICU-Metaanalysis.jags", data = data,
                         n.chains = 3, inits = inits)

ndraws <- 5000
burnin <- 10000
update(jags_model, burnin)
jags_output <- coda.samples(jags_model,
                            variable.names = c('twentypct', 'eightypct', 'sigmag', 'ORf'),
                            n.iter=ndraws)
jags_df <- ggs(jags_output)

### Check convergence
ggs_traceplot(jags_df) + facet_wrap(~Parameter, scales = "free_y")

jags_df %>%
    group_by(Parameter) %>%
    summarize(mean = mean(value),
              sd = sd(value),
              `2.5%` = quantile(value, .025),
              median = median(value),
              `97.5%` = quantile(value, .975))

### Now stan

stan_fit <- stan(file = 'Ch14-ICU-Metaanalysis.stan', data = data, chains = 3, 
				iter = (ndraws+burnin), warmup = burnin,
                init = rep(list(inits), 3))
                
stan_df <- ggs(stan_fit) %>%
    filter(Parameter %in% c('twentypct', 'eightypct', 'sigmag', 'ORf'))

### Check convergence
ggs_traceplot(stan_df) + facet_wrap(~Parameter, scales = "free_y")
    
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



