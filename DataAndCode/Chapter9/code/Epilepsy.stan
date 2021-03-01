data {
    int<lower=0> Nobs;
    int<lower=0> Npred;
    int<lower=0> y[Nobs];
    vector[Nobs + Npred] Trt;
    vector[Nobs + Npred] Base;
    vector[Nobs + Npred] Age;
}
transformed data {
    vector[Nobs + Npred] slogBase;
    vector[Nobs + Npred] slogAge;
    slogBase = (log(Base) - mean(log(Base[1:Nobs]))) / sd(log(Base[1:Nobs]));
    slogAge = (log(Age) - mean(log(Age[1:Nobs]))) / sd(log(Age[1:Nobs]));
}
parameters {
    vector[5] b;
}
transformed parameters {
    vector<lower=0>[Nobs + Npred] lam;
    for(j in 1:(Nobs + Npred)) {
        lam[j] = exp(b[1] + b[2]*slogBase[j] + b[3]*Trt[j] + b[4]*slogBase[j]*Trt[j] + b[5]*slogAge[j]);
    }
}

model {
    // Likelihood
    for(j in 1:Nobs) {
        y[j] ~ poisson(lam[j]);
    }
    // priors on b
    for(i in 1:5) {
        b[i] ~ normal(0,1);
    }
}
generated quantities {
    // Predictions for new covariate combinations
    vector[Npred] ypred;
    for(j in (Nobs + 1):(Nobs + Npred)) {
        ypred[j - Nobs] = poisson_rng(lam[j]);
    }
}