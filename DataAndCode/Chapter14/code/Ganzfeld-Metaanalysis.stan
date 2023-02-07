data {
    int nstudy;
    int n[nstudy];
    int y[nstudy];
    real theta0;
    real u;
    real thetamax;
}
transformed data{
    real mu0;

    mu0 = log(theta0)/(1+log(theta0));
}
parameters {
    real mu;
//    real z;
    vector[nstudy] w;
    real<lower=0> sigt;
}
transformed parameters {
    real c;
    real ten;
    real med;
    real ninety;
    real sigmax;
    vector[nstudy] theta;
    
    c = (logit(u) - logit(theta0))/1.645;
    for(i in 1:nstudy) {
      theta[i] =  inv_logit(w[i]);
    }
    sigmax = ((logit(thetamax) - logit(theta0))/1.28); // + z;
    ten = inv_logit(mu - (1.28*sigt));
    med = inv_logit(mu);
    ninety =  inv_logit(mu + (1.28*sigt));
}
model {
   for (i in 1:nstudy) { 
      w[i] ~ normal(mu, sigt);
      y[i] ~ binomial(n[i], theta[i]);
   }
   mu ~ normal(mu0, c); 
   sigt ~ uniform(0, sigmax);
//   z ~ uniform(0, 0.0001);
}
generated quantities {
    vector[6] prob;

    prob[1] = med > 0.25;
    prob[2] = med > 0.28;
    prob[3] = med > 0.30;
    prob[4] = ten > 0.25;
    prob[5] = ten > 0.28;
    prob[6] = ten > 0.30;
}

