library(rstan)
library(coda)
library(ggmcmc)


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
inits = list(beta=c(0), # lamstar = 0.05,
             lam=rep(.01, 23))



### Now stan
stan_fit <- stan(file = 'PHRegression.stan', data = data, chains = 3,
                 init = rep(list(inits), 3))
stan_df <- ggs(stan_fit) %>% 
    filter(Parameter %in% c('HR', 'beta'))
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


