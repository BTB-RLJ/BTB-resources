data {
    int<lower=0> y[2];
    vector[2] M;
    vector[2] Sex;
}
parameters {
    vector[2] b;
}
transformed parameters {
    vector<lower=0>[2] llam;
    vector<lower=0>[2] lam;
    real<lower=0> RR;
    for(i in 1:2) {
        llam[i] = exp(log(M[i]) + b[1] + b[2]*Sex[i]);
    }
    for(i in 1:2) {
        lam[i] = exp(b[1] + b[2]*Sex[i]);
    }
    RR = lam[1]/lam[2];
}
model {
    for(i in 1:2) {
        y[i] ~ poisson(llam[i]);
    }
    for(i in 1:2) {
        b[i] ~ normal(0,1);
    }
}