data {
    int B[37];
    vector[37] y;
}
parameters {
    real<lower=0, upper=100> mu;
    vector[6] gamma;
    real<lower=0, upper=50> sigma;
    real<lower=0, upper=100> sigmag;
}
model {
    for(i in 1:37) {
        y[i] ~ normal(gamma[B[i]], sigma);
    }
    gamma ~ normal(mu, sigmag);
}
generated quantities {
   int prob;
   vector[2] R;

   R[1] = sigmag/sigma;
   R[2] = gamma[5]/gamma[4];
   prob = (gamma[5] > gamma[4]);  // returns 1 if gamma[5]>gamma[4]
}
