data {
    int N;
    int nID;
    int MyID[N];
    vector[N] time;
    vector[N] trt;
    int<lower=0, upper=1> y[N];
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
    vector[nID] gamma;
    real<lower=0, upper=10> gamma_sd;
    vector[4] beta;
}
transformed parameters {
    vector<lower=0,upper=1>[N] p;
    for(i in 1:N) {
       p[i] = inv_logit(gamma[MyID[i]] + beta[2]*(trt[i] * z[i]) + beta[3]*sTime[i] + beta[4]*(trt[i]*sTime[i]*z[i]));
    }
}
model {
    for (i in 1:N) {
        y[i] ~ bernoulli(p[i]);
    }
    for (j in 1:nID) {
      gamma[j] ~ normal(beta[1], gamma_sd);
    }
   for (k in 1:4) {
      beta[k] ~ normal(0, 1);
   }
}
