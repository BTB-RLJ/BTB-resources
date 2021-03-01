data {
    int<lower=0> y[4];
    int<lower=0> n[4];
    vector[4] sDose;
}
parameters {
    vector[2] beta;
}
transformed parameters {
    real p[4];
    for(i in 1:4){
         p[i] = inv_logit(beta[1] + beta[2] * sDose[i] );
    }
}
model {
   for(i in 1:4) {
      y[i] ~ binomial(n[i], p[i]);
   }
    // Prior
   for(i in 1:2) {
      beta[i] ~ uniform(-10,10);
   }
}