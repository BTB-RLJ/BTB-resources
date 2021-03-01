data {
    vector[10] A;
    int S[10];
    int<lower=0> y[10];
    vector<lower=0>[10] M;
    real b10;
    real<lower=0> tau10;
    real b20;
    real<lower=0> tau20;
}
transformed data {
    vector[10] sA;
    vector[10] sAsq;
    sA = (A - mean(A)) / sd(A);
    for(i in 1:10) {
        sAsq[i] = sA[i] ^ 2;
    }
}
parameters {
    vector[4] beta;
    vector[2] lam;
}
transformed parameters {
    vector[10] m;
    vector[5] RR;
    vector[10] lambda;
    vector[10] llambda;
    vector[2] b;
    b[1] = log(lam[1]);
    b[2] = log(lam[2]) - b[1];
    for(i in 1:10) {
         m[i] = b[1]+b[2]*S[i]+beta[1]*sA[i]+beta[2]*sAsq[i] + beta[3]*sA[i]*S[i] +beta[4]*sAsq[i]*S[i];
         lambda[i] = exp(m[i]);
         llambda[i] = exp(log(M[i]) + log(lambda[i]));
    }
    for(i in 1:5) {
        RR[i] = lambda[i] /lambda[i+5];
    }

}

model {
   for(i in 1:10) {
       y[i] ~ poisson(llambda[i]);
   }
   lam[1] ~ lognormal(b10, 1 / sqrt(tau10));
   lam[2] ~ lognormal(b20, 1 / sqrt(tau20));
   for(i in 1:4) {
       beta[i] ~ uniform(-10, 10);
   }
}