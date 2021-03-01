data {
    int<lower=0> n[3];
    vector[n[1]] low;
    vector[n[2]] normal;
    vector[n[3]] high;
}
parameters {
    vector[3] mu;
    real<lower=0> tau;
}
transformed parameters {
    real<lower=0>sigma;
    sigma = 1/sqrt(tau);
}
model {
    // Likelihood
    for(i in 1:n[1]) {
        low[i] ~ lognormal(mu[1],sigma);
    }
    for(i in 1:n[2]) {
        normal[i] ~ lognormal(mu[2],sigma);
    }
    for(i in 1:n[3]) {
        high[i] ~ lognormal(mu[3],sigma);
    }
    // Prior
    mu[1]~normal(4.87, 1 / sqrt(347.12));
    mu[2]~normal(5.39, 1 / sqrt(357.42));
    mu[3]~normal(6.40, 1 / sqrt(166.67));
    tau~gamma(0.001, 0.001);
}
generated quantities {
    // Note that not every generated quantity in BUGS code is included here
    vector<lower=0>[3] med;
    for(i in 1:3) {
        med[i] = exp(mu[i]);
    }
}