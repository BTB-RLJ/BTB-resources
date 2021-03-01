data {
    int<lower=0> N;
    vector[N] DBP;
    vector[N] Bdp;
    vector[N] Wt;
    vector[N] Ht;
    vector[N] Trt;
}
parameters {
    vector[5] b;
    real gamma;
}
transformed parameters {
    real<lower=0>sigma;
    sigma = 1/sqrt(exp(gamma));
}
model {
    // Likelihood
    for(i in 1:N) {
        DBP[i] ~ normal(b[1] + b[2]*(Trt[i] -1) + b[3]* (Bdp[i] - 72.4) + b[4] * (Ht[i] - 68.9) + b[5] * (Wt[i] - 205.9), sigma);
    }
    // No need for priors since they are flat
}
generated quantities {
    vector[4] TrtEffect;
    vector[4] M;
    TrtEffect[1] = b[2];
    TrtEffect[2] =10*b[3];
    TrtEffect[3] = 6*b[4];
    TrtEffect[4] = 40*b[5];
    M[1] = b[1] + b[2] + (80-72.4)*b[3]; 
    M[2] = b[1] + b[2] + (70- 72.4)*b[3] ;
    M[3] = b[1] + (80-72.4)*b[3];
    M[4] = b[1] +  (70- 72.4)*b[3]; 
}