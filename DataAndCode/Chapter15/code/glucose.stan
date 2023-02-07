data {
    int n[2];
    vector[n[1] + n[2]] y;
}
parameters {
    vector<lower=0>[2] tau;
    vector[2] mu;
}
transformed parameters {
    vector<lower=0>[2] sigma;
    for(i in 1:2) {
       sigma[i] = 1 / sqrt(tau[i]); 
    }
    
}
model {
    for(i in 1:n[1]) {
        y[i] ~ normal(mu[1], sigma[1]);
    }
    for (i in (n[1] +1):(n[1] + n[2])) {
        y[i] ~ normal(mu[2], sigma[2]);
    }
   mu[1] ~ normal( 100, 1 / sqrt(0.001));
   mu[2] ~ normal( 250, 1 / sqrt(0.0001));
   tau[1] ~ gamma(.1,.1);
   tau[2] ~ gamma(.1,.1);
}
generated quantities {
    real AUC;
    AUC = Phi((mu[2] - mu[1])/sqrt(sigma[1] ^ 2 + sigma[2] ^ 2));
}
