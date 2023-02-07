library(rstan)
library(rjags)
library(coda)
library(ggmcmc)

options(mc.cores = parallel::detectCores()) # Allows Stan to use multiple cores

DentalData <- read.table('../data/DentalData.txt', header=TRUE)

# Data
data <- list(nObs = nrow(DentalData),
			nID = length(unique(DentalData$ID1)),
			dist = DentalData$Distance,
            Age = DentalData$Age,
            ID = DentalData$ID1,
            Male = DentalData$Male)

inits = function(){
	list(beta=rnorm(3, 0, 1), sig=runif(1, 0.5, 2.0), sigg=runif(2,0.5,2.0), rhog = runif(1, 0.0, 0.5))
}

jags_model <- jags.model(file = "dental.jags", data = data, inits=inits,
                         n.chains = 3)

ndraws <- 5000
burnin <- 5000
update(jags_model, burnin)
jags_output <- coda.samples(jags_model,
                            variable.names = c('beta', 'sig', 
                                               'sigg', 'rhog'),
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
stan_fit <- stan(file = 'dental.stan', data = data, chains = 3, 
				iter = 10000, warmup = 5000,
                init = inits)

stan_df <- ggs(stan_fit) %>%
    filter(Parameter %in% c('beta[1]', 'beta[2]', 'beta[3]',
                            'sig', 'sigg[1]', 'sigg[2]', 'rhog'))

### Check convergence
ggs_traceplot(stan_df) + facet_wrap(~Parameter, scales = "free_y")

### Can check sampling and divergent transitions
# pairs(stan_fit, pars=c('beta[1]', 'beta[2]', 'beta[3]',
                            'sig', 'sigg[1]', 'sigg[2]', 'rhog'))

stan_df %>%
    group_by(Parameter) %>%
    summarize(mean = mean(value),
              sd = sd(value),
              `2.5%` = quantile(value, .025),
              median = median(value),
              `97.5%` = quantile(value, .975))

Results <- bind_rows(mutate(jags_df, sampler = "JAGS"),
                mutate(stan_df, sampler = "STAN"))
ggplot(Results, aes(x = value, color = sampler)) +
    geom_line(stat = "density") +
    facet_wrap(~Parameter, scales = "free")
