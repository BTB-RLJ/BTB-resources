data {
    int<lower=0> n[2];
    int<lower=0> Le[2];
    int<lower=0> High[2];
}
parameters {
    vector[2] theta;
}
transformed parameters {
    vector[2] beta;
    beta[1] = logit(theta[2]);
    beta[2] = logit(theta[1]) - beta[1];
}
model {
   for (i in 1:2)  {
      Le[i] ~ binomial (n[i], inv_logit(beta[1] + beta[2] * High[i]));
   }
    // Prior
    for (i in 1:2) {
      theta[i] ~ uniform(0,1);
   }
}