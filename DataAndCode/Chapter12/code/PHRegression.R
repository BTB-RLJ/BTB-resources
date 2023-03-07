library(rstan)
library(rjags)
library(coda)
library(ggmcmc)

options(mc.cores = parallel::detectCores())

# Data
data <- list(n=33, k=23, w=0.01,
             y=c(65,156,100,134,16,108,121,4,39,143,56,26,22,1,1,5,65,
                 56,65,17,7,16,22,3,4,2,3,8,4,3,30,4,43),	
             AGgroup=c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                       0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
             delta=c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1),
             d=c(0,1,2,3,4,5,7,8,16,17,22,26,30,39,43,
                 56,65,100,108,121,134,143,156))
                 
inits = function(){
	list(beta=rnorm(1, 0, 1), lam=rgamma(23, 0.01, 1))
}

jags_model <- jags.model(file = "PHRegression.txt", data = data,
                         n.chains = 3, inits = inits)

ndraws <- 3000
burnin <- 5000
update(jags_model, burnin)
jags_output <- coda.samples(jags_model,
                            variable.names = c('HR', 'beta'),
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

ggplot(jags_df, aes(x = value)) +
    geom_line(stat = "density") +
    facet_wrap(~Parameter, scales = "free")


### Now stan
stan_fit <- stan(file = 'PHRegression.stan', data = data, chains = 3, iter=8000, warmup=5000,
                 init = inits)
                 
stan_df <- ggs(stan_fit) %>% 
    filter(Parameter %in% c('HR', 'beta'))
    
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

