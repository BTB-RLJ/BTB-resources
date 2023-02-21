data {
    int m;
    vector[m] LOR;
    vector<lower=0>[m] se;
}
parameters {
    vector[m] theta;
    real<lower=-2, upper=2> mu;
    real<lower=0, upper=2> sigmag;
}
model {
    for(i in 1:m) {
        LOR[i] ~ normal(theta[i], se[i]);
    }
    for(i in 1:m) {
        theta[i] ~ normal(mu, sigmag);
    }
//    mu ~ normal(0, sqrt(1000));
//    sigmag ~ uniform(0, 10);
}
generated quantities {
   real medOR;
   int prob;
   real eightypct;
   real twentypct;
   real muf;
   real ORf;
   vector[m] OR;
   medOR = exp(mu);
   prob = medOR > 1;
// 20th and 80th percentiles of distribution of ORs
   eightypct = exp(mu + 0.8416*sigmag);
   twentypct = exp(mu - 0.8416*sigmag);
// Predictive density of ORs
   muf = normal_rng(mu, sigmag); 
   ORf = exp(muf);
   for (i in 1:m) { 
      OR[i] = exp(theta[i]);
   }
}
