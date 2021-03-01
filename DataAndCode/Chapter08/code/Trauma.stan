data {
    int<lower=0> n;
    int<lower=0> death[n];
    vector[n] ISS;
    vector[n] TI;
    vector[n] RTS;
    vector[n] AGE;
    vector<lower=0>[6] a;
    vector<lower=0>[6] b;
    matrix[6,6] Xtinv;
}
parameters {
    real<lower=0,upper=1> tildetheta[6];
}
transformed parameters {
    vector[6] beta;
    vector[6] v;
    vector[n] theta;
    for(i in 1:6) {
      v[i] = log(tildetheta[i]/(1-tildetheta[i]));
   }
   for(i in 1:6) {
       beta[i] = dot_product(Xtinv[i,1:6], v);
   }
   for(i in 1:n) {
      theta[i] = inv_logit(beta[1] + beta[2]*ISS[i] + beta[3]*RTS[i] + beta[4]*AGE[i] + beta[5]*TI[i] + beta[6]*AGE[i]*TI[i]);
   }
   
}
model {
   for(i in 1:n) {
      death[i] ~ bernoulli(theta[i]);
   }
    // Prior
   for(i in 1:6) {
      tildetheta[i] ~ beta(a[i],b[i]);
   }
}