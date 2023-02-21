> ndraws <- 4000
> burnin <- 6000
> update(jags_model, burnin)
  |**************************************************| 100%
> jags_output <- coda.samples(jags_model,
+                             variable.names = c('beta', 'gamma_sd'),
+                             n.iter=ndraws)
  |**************************************************| 100%
> jags_df <- ggs(jags_output)
> jags_df %>%
+     group_by(Parameter) %>%
+     summarize(`2.5%` = quantile(value, .025),
+               median = median(value),
+               `97.5%` = quantile(value, .975))
# A tibble: 5 × 6
  Parameter   mean    sd `2.5%` median `97.5%`
  <fct>      <dbl> <dbl>  <dbl>  <dbl>   <dbl>
1 beta[1]   -3.10  0.324  -3.75 -3.09  -2.49  
2 beta[2]   -0.895 0.343  -1.58 -0.885 -0.232 
3 beta[3]   -1.60  0.175  -1.95 -1.59  -1.27  
4 beta[4]   -0.521 0.286  -1.08 -0.524  0.0364
5 gamma_sd   3.85  0.355   3.20  3.85   4.57  

> stan_fit <- stan(file = 'Toenail.stan', data = data, chains = 3, iter = 10000, warmup = 6000,
+                  init = rep(list(inits), 3), 
+                  pars=c("gamma", "p"), include=FALSE) # Do not save "gamma" or "p" vectors
> stan_df = ggs(stan_fit)
> stan_df %>% 
+     group_by(Parameter) %>%
+     summarize(mean = mean(value),
+               sd = sd(value),
+               `2.5%` = quantile(value, .025),
+               median = median(value),
+               `97.5%` = quantile(value, .975))
# A tibble: 5 × 6
  Parameter   mean    sd `2.5%` median `97.5%`
  <fct>      <dbl> <dbl>  <dbl>  <dbl>   <dbl>
1 beta[1]   -3.10  0.324  -3.76 -3.09  -2.50  
2 beta[2]   -0.920 0.334  -1.57 -0.918 -0.267 
3 beta[3]   -1.60  0.171  -1.95 -1.60  -1.28  
4 beta[4]   -0.522 0.281  -1.08 -0.516  0.0217
5 gamma_sd   3.89  0.351   3.25  3.87   4.63  

