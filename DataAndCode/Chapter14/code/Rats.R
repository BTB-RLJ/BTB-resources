library(rstan)
library(rjags)
library(coda)
library(ggmcmc)

# Data
data <- list(N = 30,
             T = 5,
             y = structure(c(151, 145, 147, 155, 135, 159, 141, 159, 177, 134, 
                             160, 143, 154, 171, 163, 160, 142, 156, 157, 152, 154, 139, 146, 
                             157, 132, 160, 169, 157, 137, 153, 199, 199, 214, 200, 188, 210, 
                             189, 201, 236, 182, 208, 188, 200, 221, 216, 207, 187, 203, 212, 
                             203, 205, 190, 191, 211, 185, 207, 216, 205, 180, 200, 246, 249, 
                             263, 237, 230, 252, 231, 248, 285, 220, 261, 220, 244, 270, 242, 
                             248, 234, 243, 259, 246, 253, 225, 229, 250, 237, 257, 261, 248, 
                             219, 244, 283, 293, 312, 272, 280, 298, 275, 297, 350, 260, 313, 
                             273, 289, 326, 281, 288, 280, 283, 307, 286, 298, 267, 272, 285, 
                             286, 303, 295, 289, 258, 286, 320, 354, 328, 297, 323, 331, 305, 
                             338, 376, 296, 352, 314, 325, 358, 312, 324, 316, 317, 336, 321, 
                             334, 302, 302, 323, 331, 345, 333, 316, 291, 324), .Dim = c(30, 
                                                                                         5)),
             x = c(8.0, 15.0, 22.0, 29.0, 36.0),
             xbar = 22)
inits = list(alpha=rnorm(30, 200, 50), beta=rnorm(30, 0, 3), 
			mu_alpha=rnorm(1, 200, 100), mu_beta=rnorm(1, 0, 6), 
			alpha_tau=rgamma(1, 1, 1), beta_tau=rgamma(1, 1, 1),
			tau.c=rgamma(1, 1, 1))
			
jags_model <- jags.model(file = "Rats.jags", data = data,
                         n.chains = 3, inits = inits)

ndraws <- 2000
burnin <- 2000
update(jags_model, burnin)
jags_output <- coda.samples(jags_model,
                            variable.names = c('mu_alpha', 'mu_beta'),
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

stan_fit <- stan(file = 'Rats.stan', data = data, chains = 3,
                 init =rep(list(inits), 3))
stan_df <- ggs(stan_fit) %>%
    filter(Parameter %in% c('mu_alpha', 'mu_beta'))
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


