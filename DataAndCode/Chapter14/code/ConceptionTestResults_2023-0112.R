> inits <- list(mu=runif(1,30,70), sigma=runif(1, 0.5, 1.5), sigmag=runif(1, 1.0, 2.0))
> 
> jags_model <- jags.model(file = "Ch14-Conception.jags", data = data,
+                          n.chains = 3, inits = inits)
  |++++++++++++++++++++++++++++++++++++++++++++++++++| 100%
> 
> ndraws <- 10000
> burnin <- 5000
> update(jags_model, burnin)
  |**************************************************| 100%
> jags_output <- coda.samples(jags_model,
+                             variable.names = c('mu', 'gamma', 'sigma', 'sigmag', 'R', 'prob'),
+                             n.iter=ndraws)
  |**************************************************| 100%
> jags_df <- ggs(jags_output)
> jags_df %>%
+     group_by(Parameter) %>%
+     summarize(mean = mean(value),
+               sd = sd(value),
+               `2.5%` = quantile(value, .025),
+               median = median(value),
+               `97.5%` = quantile(value, .975))
# A tibble: 12 × 6
   Parameter   mean    sd  `2.5%` median `97.5%`
   <fct>      <dbl> <dbl>   <dbl>  <dbl>   <dbl>
 1 gamma[1]  46.3   6.92  32.2    46.5     58.9 
 2 gamma[2]  58.6   8.49  43.1    57.9     76.9 
 3 gamma[3]  55.5   5.34  45.1    55.4     66.5 
 4 gamma[4]  45.3   7.14  30.9    45.5     58.3 
 5 gamma[5]  62.6   6.23  51.1    62.5     75.0 
 6 gamma[6]  55.2   4.47  46.5    55.2     64.2 
 7 mu        53.9   6.44  40.8    53.9     67.1 
 8 prob       0.958 0.200  0       1        1   
 9 R[1]       0.739 0.523  0.0574  0.635    2.06
10 R[2]       1.43  0.320  0.986   1.38     2.16
11 sigma     16.6   2.21  12.9    16.4     21.5 
12 sigmag    11.9   8.11   1.04   10.3     32.1 

Warning messages:
1: There were 37 divergent transitions after warmup. See
https://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
to find out why this is a problem and how to eliminate them. 
2: Examine the pairs() plot to diagnose sampling problems
 
>                 
> stan_df <- ggs(stan_fit) 
> stan_df %>%
+     group_by(Parameter) %>%
+     summarize(mean = mean(value),
+               sd = sd(value),
+               `2.5%` = quantile(value, .025),
+               median = median(value),
+               `97.5%` = quantile(value, .975))
# A tibble: 12 × 6
   Parameter   mean    sd `2.5%` median `97.5%`
   <fct>      <dbl> <dbl>  <dbl>  <dbl>   <dbl>
 1 gamma[1]  46.1   6.77  32.5   46.2     58.6 
 2 gamma[2]  58.8   8.53  43.3   58.1     77.3 
 3 gamma[3]  55.5   5.31  45.2   55.4     66.5 
 4 gamma[4]  45.1   7.00  30.9   45.3     57.9 
 5 gamma[5]  62.9   6.27  51.1   62.8     75.4 
 6 gamma[6]  55.2   4.53  46.4   55.2     64.4 
 7 mu        54.0   6.60  40.2   53.9     67.4 
 8 prob       0.969 0.174  0      1        1   
 9 R[1]       0.757 0.520  0.115  0.651    2.12
10 R[2]       1.44  0.312  0.992  1.39     2.17
11 sigma     16.5   2.20  12.8   16.3     21.4 
12 sigmag    12.2   8.05   2.05  10.6     33.3 
> stan_fit <- stan(file = 'Ch14-Conception.stan', data = data, chains = 3, 
+ 				iter = 25000, warmup = 20000, 
+ 				control = list(adapt_delta=0.99),
+                 init = rep(list(inits), 3))

