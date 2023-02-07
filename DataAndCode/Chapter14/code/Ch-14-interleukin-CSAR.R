library(rstan)
library(rjags)
library(coda)
library(ggmcmc)

options(mc.cores = parallel::detectCores()) # Allows Stan to use multiple cores

ILdata = read.table('../data/IL-1betaData.txt', header=TRUE)

u.ID = unique(ILdata$ID)
n.ID = length(u.ID)

#	Convert data to matrix format for autoregression program
ILmatrix = matrix(NA, nrow=n.ID, ncol=6)
#	Put data in matrix format for Wes's program
for (i in 1:n.ID) {
	ILmatrix[i,] = ILdata[ILdata[,1]==u.ID[i],2]
}
             
inits = function() {
	list(sig=runif(1, 0.5, 1.5), sigg=runif(1, 0.5, 1.5), sigw=runif(1, 0.5, 1.5), rhow=runif(1, 0.25, 0.75))
}

data <- list(N=n.ID, logy = log(ILmatrix)) # Fit logarithms of IL-1beta values
jags_model <- jags.model(file = "Ch-14-interleukin-CSAR.jags", data = data, n.chains = 3)

ndraws <- 10000
burnin <- 6000
update(jags_model, burnin)
jags_output <- coda.samples(jags_model,
                            variable.names = c('med', 'ten', 'ninety', 
                            					'sig', 'sigg', 'sigw', 'rho1', 'rho2', 'rhow'),
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
stan_fit <- stan(file = 'Ch14-interleukin-CSAR.stan', data = data, chains = 3, iter=10000, warmup=6000,
                                  init = inits)

stan_df <- ggs(stan_fit) %>%
	filter(Parameter %in% c('med', 'ten[1]', 'ten[2]', 'ninety[1]', 'ninety[2]', 
                            'sig', 'sigg', 'sigw', 'rho1', 'rho2', 'rhow'))

### Check poor sampling and divergent transitions
pairs(stan_fit, pars=c("beta1", "rhow", "sig", "sigg", "sigw"))

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
