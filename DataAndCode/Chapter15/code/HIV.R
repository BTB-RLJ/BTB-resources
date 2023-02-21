library(rstan)
library(rjags)
library(coda)
library(ggmcmc)

options(mc.cores = parallel::detectCores()) # Allows Stan to use multiple cores

data <- list(y = c(64, 309,3122183 ), n= 3122556,
             a = 15, b = 94092,
             a1 = 142, b1 = 2,
             a2 = 1363, b2 = 3)
             
inits <- function(){
	list(pi=runif(1, 0.1, 0.6), se=runif(1, 0.4, 0.9), sp=runif(1, 0.4, 0.9))
}
             
jags_model <- jags.model(file = "HIV.jags", data = data, inits = inits,
                         n.chains = 3)

ndraws <- 8000
burnin <- 4000

update(jags_model, burnin)

jags_output <- coda.samples(jags_model,
                            variable.names = c('se', 'sp', 'pi', 'ppv', ' pdgneg'),
                            n.iter=ndraws)

jags_df <- ggs(jags_output)

### Can check convergence
ggs_traceplot(jags_df) + facet_wrap(~Parameter, scales = "free_y")

jags_df %>%
    group_by(Parameter) %>%
    summarize(mean = mean(value),
              sd = sd(value),
              `2.5%` = quantile(value, .025),
              median = median(value),
              `97.5%` = quantile(value, .975))

### Now stan
stan_fit <- stan(file = 'HIV.stan', data = data, chains = 3,
				iter = (ndraws+burnin), warmup = burnin,
                init = inits)

stan_df <- ggs(stan_fit) %>%
    filter(Parameter %in% c('se', 'sp', 'pi', 'ppv', 'pdgneg'))

### Can check convergence
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



