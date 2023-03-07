library(rjags)
library(rstan)
library(coda)
library(ggmcmc)

options(mc.cores = parallel::detectCores())

df <- read.table("../data/Larynx-Cancer-Data.txt", header = TRUE) 
colnames(df) <- c("stage", "time", "logtime", "age", "yr", "cens_time", "log_cens_time")

# Data
stage <- df$stage
t <- df$time
age <- df$age
yr <- df$yr
c <- df$cens_time
is.censored <- is.na(t)
c[!is.censored] <- t[!is.censored] + 1

data <- list(oldn = 90, predn = 16,
             stage = stage, t = t[1:90], predt = rep(NA, 16), age = age,
             yr = yr, is.censored = is.censored[1:90], c = c[1:90])
inits = list(beta=c(0,0,0,0,0,0),tau=1)
inits = function(){
	list(beta=rnorm(6, 0, 1), tau=runif(1, 0.5, 2.0))
}

jags_model <- jags.model(file = "LarynxLogNormal.jags", data = data,
                         n.chains = 3, inits = inits)

ndraws <- 8000
burnin <- 5000
update(jags_model, burnin)
jags_output <- coda.samples(jags_model,
                            variable.names = c('beta'),
                            n.iter=ndraws)
jags_df <- ggs(jags_output)

### Show convergence
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
df <- df[1:90,]
stage <- df$stage
t <- df$time
age <- df$age
yr <- df$yr
c <- df$cens_time
is.censored <- is.na(t)
stan_data <- list(Nobs = sum(!is.censored),
                  Ncens = sum(is.censored),
                  tobs = t[!is.censored],
                  tcens = c[is.censored],
                  age = age,
                  yr = yr,
                  stageobs = stage[!is.censored],
                  ageobs = age[!is.censored],
                  yrobs = yr[!is.censored],
                  stagecens = stage[is.censored],
                  agecens = age[is.censored],
                  yrcens = yr[is.censored])

stan_fit <- stan(file = 'LarynxLogNormal.stan', data = stan_data, chains = 3, iter=8000, warmup=5000,
                 init = inits)

stan_df <- ggs(stan_fit, family =  "beta")

### Show convergence
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
