library(rstan)
library(rjags)
library(coda)
library(ggmcmc)

options(mc.cores = parallel::detectCores()) # Allows Stan to use multiple cores

CowData <- read.table("../data/Ch14-Conception.txt", header=TRUE)
data <- list(y=CowData$y, B=CowData$B)

inits <- function(){
	list(mu=runif(1,30,70), sigma=runif(1, 0.5, 1.5), sigmag=runif(1, 1.0, 2.0))
}

jags_model <- jags.model(file = "Ch14-Conception.jags", data = data,
                         n.chains = 3, inits = inits)

ndraws <- 10000
burnin <- 5000
update(jags_model, burnin)
jags_output <- coda.samples(jags_model,
                            variable.names = c('mu', 'gamma', 'sigma', 'sigmag', 'R', 'prob'),
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

stan_fit <- stan(file = 'Ch14-Conception.stan', data = data, chains = 3, 
				iter = (ndraws+burnin), warmup = burnin, 
				control = list(adapt_delta=0.99),
                init = inits)
                
stan_df <- ggs(stan_fit) 

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



