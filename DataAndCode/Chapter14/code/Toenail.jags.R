library(rjags)
library(coda)
library(ggmcmc)

df <- read.table('ToenailData.txt', header = TRUE)

# Data
data <- list(y = df$y,
             trt = df$Trt,
             time = df$Time,
             N = nrow(df),
             MyID = df$ID,
             nID = max(unique(df$ID)))
#inits = list(beta=c(0,0,0),alpha0=0,sigma.alpha=1)

inits = list(beta=rnorm(4, 0, 1), gamma=rep(0, data$nID), gamma_sd=1)
jags_model <- jags.model(file = "Toenail.txt", data = data,
                         n.chains = 2, inits = inits)

ndraws <- 4000
burnin <- 6000
update(jags_model, burnin)
jags_output <- coda.samples(jags_model,
                            variable.names = c('beta', 'gamma_sd'),
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

