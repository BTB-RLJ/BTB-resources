library(rstan)
library(rjags)
library(coda)
library(ggmcmc)


options(mc.cores = parallel::detectCores()) # Allows Stan to use multiple cores

### Read data
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
inits <- function(){
	list(mu=rnorm(3, c(4.87,5.39,6.4), 1),tau=rgamma(1,1,1))
}

jags_model <- jags.model(file = "DiasorinANOVA.jags", data = data,
                         n.chains = 3, inits = inits)

ndraws <- 2000
burnin <- 2000

update(jags_model, burnin)

jags_output <- coda.samples(jags_model,
                            variable.names = c('med', 'sigma', 'P', 'diff21', 'diff31', 'diff32', 
                            'prob21', 'prob31', 'prob32'),
                            n.iter=ndraws)

jags_df <- ggs(jags_output)

### Check convergence
ggs_traceplot(jags_df) + facet_wrap(~Parameter, scales = 'free_y')

jags_df %>%
    group_by(Parameter) %>%
    summarize(mean = mean(value),
              sd = sd(value),
              `2.5%` = quantile(value, .025),
              median = median(value),
              `97.5%` = quantile(value, .975))

### Now stan
stan_fit <- stan(file = 'DiasorinANOVA.stan', data = data, chains = 3,
					pars=c('med', 'sigma', 'P', 'diff21', 'diff31', 'diff32', 
                            'prob21', 'prob31', 'prob32'),
					init=inits, iter=4000, warmup=2000)

stan_df <- ggs(stan_fit)

### Check convergence
ggs_traceplot(stan_df) + facet_wrap(~Parameter, scales = 'free_y')

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


