library(rstan)
library(rjags)
library(coda)
library(ggmcmc)

data <- list(y=c(14,14,11,13,55,22,12,93,22,33,66,30,16,42,59,
                 
                 16,6,123,15,16,14,14,13,30,143,6,10,53,42,28,7,
                 
                 13,19,11,74,20,10,24,29,4,6,12,65,26,39,7,32,3,
                 
                 302,13,26,10,70,13,15,51,6,0,10, 
                 
                 NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),
             
             Trt = c( 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                      
                      0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                      
                      0, 0, 0, 0, 0, 0, 0, 0, 1, 1,
                      
                      1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                      
                      1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                      
                      1, 1, 1, 1, 1, 1, 1, 1, 1, 
                      
                      1,1,1,0,0,0, 1,1,1,0,0,0),
             
             Base = c( 11, 11,  6,  8, 66, 27, 12, 52, 23, 10,
                       
                       52, 33, 18, 42, 87, 50, 18,111, 18, 20,
                       
                       12,  9, 17, 28, 55,  9, 10, 47, 76, 38,
                       
                       19, 10, 19, 24, 31, 14, 11, 67, 41,  7,
                       
                       22, 13, 46, 36, 38,  7, 36, 11,151, 22,
                       
                       41, 32, 56, 24, 16, 22, 25, 13, 12,  
                       
                       22,22,22, 22,22,22, 50,50,50,50,50,50 ),
             
             Age  = c(31,30,25,36,22,29,31,42,37,28,
                      
                      36,24,23,36,26,26,28,31,32,21,
                      
                      29,21,32,25,30,40,19,22,18,32,
                      
                      20,30,18,24,30,35,27,20,22,28,
                      
                      23,40,33,21,35,25,26,25,22,32,
                      
                      25,35,21,41,32,26,21,36,37,  
                      
                      20, 28, 35, 20,28, 35, 20,28,35,20,28,35))

### Note that data is ordered by observed data and then covariates for which we want to predict
data <- c(data, Nobs = sum(!is.na(data$y)), Npred = sum(is.na(data$y)))
data$y <- data$y[!is.na(data$y)]
jags_inits <- list(b = c(0,0,0,0,0)) 

jags_model <- jags.model(file = "Epilepsy.txt", data = data,
                         n.chains = 3, inits = jags_inits)

ndraws <- 1000
burnin <- 1000
update(jags_model, burnin)
jags_output <- coda.samples(jags_model,
                            variable.names = c('b'),
                            n.iter=ndraws)
jags_df <- ggs(jags_output)
jags_df %>%
    group_by(Parameter) %>%
    summarize(mean = mean(value),
              median = median(value),
              `2.5%` = quantile(value, .025),
              `97.5%` = quantile(value, .975))

### Now stan
stan_fit <- stan(file = 'Epilepsy.stan', data = data, chains = 3,
                 init = "0")
stan_df <- ggs(stan_fit, family = "b")
stan_df %>%
    group_by(Parameter) %>%
    summarize(mean = mean(value),
              median = median(value),
              `2.5%` = quantile(value, .025),
              `97.5%` = quantile(value, .975))

df <- bind_rows(mutate(jags_df, sampler = "JAGS"),
                mutate(stan_df, sampler = "STAN"))
ggplot(df, aes(x = value, color = sampler)) +
    geom_line(stat = "density") +
    facet_wrap(~Parameter, scales = "free")



