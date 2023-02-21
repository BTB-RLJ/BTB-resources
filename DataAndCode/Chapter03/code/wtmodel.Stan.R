library(rstan)
library(ggmcmc)

options(mc.cores = parallel::detectCores())  # Allows Stan to use multiple cores for parallel processing

VetBPdata = read.csv("../data/VetBP.csv", header=T)

ok = !is.na(VetBPdata$Weight_00)
y = VetBPdata$Weight_00[ok]
n = length(y)
alpha0 <- 0.1; lambda0 <- 0.1;
mu0 <- 150; tau0 <- 0.001

data = list("n"=n, "y"=y, "mu0"=mu0, "tau0"=tau0, "alpha0"=alpha0, "lambda0"=lambda0)

inits = function() {
  list(mu = rnorm(1, mean=120, sd=5), tau = runif(1, min=0.1, max=10))
}

# Run with Stan
VetBP.stan.sim = stan(file="wtmodel.stan", data=data, chains = 3, iter=6000, warmup=4000,
                  pars=c("mu", "sigma", "cv"),
                  init = inits)

stan_df <- ggs(VetBP.stan.sim) 

### Check convergence
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
