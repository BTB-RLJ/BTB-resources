data {
    int N;
    int nID;
    int MyID[N];
    vector[N] time;
    vector[N] trt;
    int<lower=0> y[N];
}
transformed data{
    int z[N];
    vector[N] sTime;
    for(i in 1:N) {
        z[i] = time[i] > 0;
    }
    sTime = (time - mean(time)) / sd(time);
}
parameters {
    vector[N] alpha;
    real alpha0;
    real<lower=0, upper=10> alpha_sd;
    vector[3] beta;
}
transformed parameters {
    vector<lower=0,upper=1>[N] p;
    for(i in 1:N) {
       p[i] = inv_logit(alpha[MyID[i]] + beta[1]*(trt[i] * z[i]) + beta[2]*sTime[i] + beta[3]*(trt[i]*sTime[i] * z[i]));
    }
}
model {
    for (i in 1:N) {
        y[i] ~ bernoulli(p[i]);
    }
    for (j in 1:nID) {
      alpha[j] ~ normal(alpha0, alpha_sd);
    }
   for (k in 1:3) {
      beta[k] ~ normal(0, 1);
   }
   alpha0 ~ normal(0, 1);
}
