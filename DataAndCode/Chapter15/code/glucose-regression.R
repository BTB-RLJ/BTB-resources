#	Glucose regression exercise
library(rstan)
library(rjags)
library(coda)
library(ggmcmc)

options(mc.cores = parallel::detectCores()) # Allows Stan to use multiple cores

Group = c(rep(1, 198), rep(2,88))
data = list(n= c(198, 88), 
y=c(82, 112,  82,  80, 121, 108, 107, 103, 156, 109,  96, 111, 128, 173, 117,  93,  79, 116, 102, 104, 114,
 146, 101, 115, 105, 105, 142,  96,  93, 140, 120, 141, 137, 100, 102, 107, 207, 109,  87, 104, 119, 101,
 118, 107, 130, 110,  96,  95, 133,  97,  92,  89, 100,  91, 127, 122, 113, 113,  81,  76,  85,  93, 107,
  98, 114, 126,  88, 105, 133, 100, 116, 120, 107,  75, 104, 110, 124, 275,  79, 111,  85,  85,  99,  72,
 136,  86,  99,  81,  87,  78,  97,  83,  88,  92, 102, 100, 109, 108,  86,  84,  95,  99, 100,  70,  78,
  87,  76,  88,  74,  89, 109, 118, 106, 110,  74,  81, 145,  71,  85,  86, 109,  96,  74,  72,  57, 182,
  90,  90,  84,  72,  83,  84,  96, 103, 108,  90,  66,  90,  97, 131, 239,  91,  88,  77, 104,  87, 100,
  79,  77, 213,  91,  74,  86, 151,  82, 116,  88,  79,  78,  90,  88,  74, 103,  76,  91,  75,  80,  92,
  78, 105,  96,  92, 107,  88,  63,  96,  87,  82,  80,  70,  84,  95,  82,  99,  84,  92,  84,  99,  77,
  75,  65, 122,  77, 137,  80,  72,  79, 107,
  87, 132, 113, 207, 186, 272, 250, 149, 198, 299, 190, 134, 183, 345, 373, 439, 306, 222, 274, 207, 110,
 283, 244, 278, 229, 297, 345, 158, 154, 217, 291, 110, 232, 314, 272, 113, 130, 183, 133, 115, 169, 219,
 223, 153, 445, 249, 165, 407, 229, 432, 109, 175, 270, 106, 150, 330,  82, 106, 232, 189, 112, 104, 108,
  91, 189, 139,  91,  88, 282, 340, 238, 330, 374, 484, 396, 370, 173, 100, 156,  89,  80, 116, 139, 227,
 124, 215, 119, 210),
age =c(37, 47, 20, 54, 62, 50, 57, 40, 51, 50, 74, 70, 60, 70, 60, 35, 35, 65, 62, 70, 28, 82, 60, 70, 32, 27, 77, 52,
 23, 62, 60, 65, 84, 42, 25, 55, 67, 60, 64, 60, 46, 60, 51, 64, 42, 58, 40, 54, 48, 52, 41, 35, 54, 64, 65, 60,
 43, 55, 70, 20, 35, 22, 52, 48, 57, 56, 58, 37, 45, 34, 30, 43, 45, 54, 50, 45, 40, 48, 48, 62, 27, 60, 21, 65,
 80, 55, 52, 49, 23, 23, 46, 20, 48, 46, 48, 48, 21, 40, 70, 26, 52, 66, 61, 32, 40, 58, 23, 31, 64, 40, 60, 47,
 20, 43, 22, 46, 48, 23, 26, 72, 21, 46, 45, 40, 65, 63, 54, 52, 35, 48, 31, 52, 40, 43, 40, 42, 68, 43, 45, 60,
 77, 44, 35, 22, 30, 25, 58, 36, 30, 50, 52, 41, 54, 80, 25, 61, 41, 40, 27, 44, 32, 48, 30, 76, 62, 61, 70, 50,
 42, 61, 65, 51, 38, 50, 29, 50, 31, 24, 60, 56, 40, 29, 48, 54, 26, 43, 39, 52, 28, 21, 49, 55, 42, 71, 33, 57,
 43, 89, 51, 48, 49, 63, 70, 58, 50, 78, 60, 53, 75, 45, 53, 54, 52, 57, 50, 54, 65, 67, 48, 46, 37, 55, 50, 60, 53, 40,
 45, 42, 55, 60, 72, 43, 52, 60, 58, 67, 63, 42, 53, 42, 54, 52, 40, 70, 47, 55, 55, 58, 51, 58, 39, 47, 60, 58,
 47, 38, 65, 50, 45, 63, 54, 53, 66, 44, 70, 52, 63, 58, 42, 28, 53, 42, 54, 64, 52, 50, 50, 60, 65, 55, 27, 63,
 47, 66, 53, 46) )

inits <- function(){
	list(b0=c(rnorm(1, 100, 10),rnorm(1, 0, 10)), b1=c(rnorm(1, 250, 10),rnorm(1, 0, 10)),
		tau = rgamma(2, 1, 1))
}

jags_model <- jags.model(file = "glucose-regression.jags", data = data, 
	inits = inits, n.chains = 3)

ndraws <- 8000
burnin <- 5000
update(jags_model, burnin)
jags_output <- coda.samples(jags_model,
                            variable.names = c('b0', 'b1', 'sigma', 'prob0', 'prob1', 
                            			'AUC', 'diffAUC', 'probAUC', 'a', 'b', 'c'),
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
stan_fit <- stan(file = 'glucose-regression.stan', data = data, chains = 3,
				pars = c('b0', 'b1', 'sigma', 'prob0', 'prob1', 
                            			'AUC', 'diffAUC', 'probAUC', 'a', 'b', 'c'),
				iter = 8000, warmup = 5000,
                init = inits)
                 
stan_df <- ggs(stan_fit) 
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


