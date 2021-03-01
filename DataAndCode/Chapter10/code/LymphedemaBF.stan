data {
    int<lower=0> n[2];
    int<lower=0> y[2];
    int<lower=0> High[2];
    vector<lower=0>[2] a;
    vector<lower=0>[2] b;
    matrix[2,2] Xtinv;
}
parameters {
    vector[2] tildetheta;
}
transformed parameters {
    vector[2] beta;
    vector[2] theta;
    beta[1] = Xtinv[1,1]*logit(tildetheta[1])+Xtinv[1,2]*logit(tildetheta[2]);
    beta[2] = Xtinv[2,1]*logit(tildetheta[1])+Xtinv[2,2]*logit(tildetheta[2]);
    for(i in 1:2) {
        theta[i] = inv_logit(beta[1] + beta[2]*High[i]);
    }
    
}
model {
   for (i in 1:2)  {
      y[i] ~ binomial (n[i], inv_logit(beta[1] + beta[2] * High[i]));
   }
    // Prior
    for (i in 1:2) {
      tildetheta[i] ~ beta(a[i],b[i]);
   }
}
generated quantities {
    vector[2] thetastar;
    vector[2] v;
    vector[2] tildethetastar;
    vector[2] tttildetheta;
    vector[2] u;
    real w;
    for(i in 1:2) {
        thetastar[i] = inv_cloglog(beta[1]+beta[2]*High[i]);
        v[i] = y[i]*(log(thetastar[i])-log(theta[i])) + (n[i]-y[i])*(log(1-thetastar[i])-log(1-theta[i]));
    }
     for (i in 1:2) {
      tildethetastar[i] = inv_cloglog(beta[1] + beta[2]*High[i]);
      tttildetheta[i] = inv_logit(beta[1] + beta[2]*High[i]);
      u[i] = (a[i]-1)*log(tildethetastar[i]) - a[i]*log(tttildetheta[i]) + b[i]*(log(1-tildethetastar[i]) - log(1-tttildetheta[i])) + beta[1]+ beta[2]*High[i];
   }
   w = exp(sum(v) + sum(u));
}