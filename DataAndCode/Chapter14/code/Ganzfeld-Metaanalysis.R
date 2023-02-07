library(rstan)
library(rjags)
library(coda)
library(ggmcmc)

options(mc.cores = parallel::detectCores())

Ganzfeld = read.table('../data/Ch14-Ganzfeld-MetaanalysisData.csv', sep=',', header=TRUE)

# Fit logarithms of IL-1beta values. 
data <- list(nstudy=nrow(Ganzfeld), y=Ganzfeld$y, n=Ganzfeld$n, theta0=0.25, u=0.8, thetamax=0.99) 

inits = function(){
	list(mu=rnorm(1, -1, 1))
}

# jags_model <- jags.model(file = "Ch14-Ganzfeld-Metaanalysis.txt", data = data, n.chains = 3)
jags_model <- jags.model(file = "Ganzfeld-Metaanalysis.txt", data = data, n.chains = 3, inits = inits)

ndraws <- 8000
burnin <- 5000
update(jags_model, burnin)
jags_output <- coda.samples(jags_model,
                            variable.names = c('med', 'ten', 'ninety'),
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
stan_fit <- stan(file = 'Ch14-Ganzfeld-Metaanalysis.stan', data = data, chains = 3, iter=8000, warmup=5000,
                 init = inits)
                 
stan_df <- ggs(stan_fit) %>% 
    filter(Parameter %in% c('med', 'ten', 'ninety'))

### Check convergence
ggs_traceplot(stan_df) + facet_wrap(~Parameter, scales = "free_y")
    
stan_df %>%
    group_by(Parameter) %>%
    summarize(mean = mean(value),
              sd = sd(value),
              `2.5%` = quantile(value, .025),
              median = median(value),
              `97.5%` = quantile(value, .975))

ggplot(stan_df, aes(x = value)) +
    geom_line(stat = "density") +
    facet_wrap(~Parameter, scales = "free")

### Compare results for JAGS and Stan
df <- bind_rows(mutate(jags_df, sampler = "JAGS"),
                mutate(stan_df, sampler = "STAN"))
ggplot(df, aes(x = value, color = sampler)) +
    geom_line(stat = "density") +
    facet_wrap(~Parameter, scales = "free")

