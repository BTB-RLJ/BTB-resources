data {
    int n[21];
    int r[21];
    vector[21] l;
    vector[21] u;
}
parameters {
    real mu;
    real<lower=.001> tau;
}
transformed parameters {
    real sigma;
    vector<lower=0,upper=1>[21] p;
    sigma = 1 / sqrt(tau);
    for(i in 1:21) {
        if(l[i] == 0)
            p[i] = Phi((log(u[i]) - mu)/sigma);
        else
            p[i] = Phi((log(u[i]) - mu)/sigma) - Phi((log(l[i]) - mu)/sigma);
    }
}
model {
    for(i in 1:21) {
       n[i] ~ binomial(r[i], p[i]); 
    }
    mu ~ normal(0, 1 / sqrt(0.001));
    tau ~ gamma(0.001,0.001)T[0.001,];
}
generated quantities {
    real med;
    real ninetypct;
    real surv;
    med = exp(mu);
    ninetypct = med*exp(1.28*sigma);
    surv = 1-Phi((log(16) - mu)/sigma);
}