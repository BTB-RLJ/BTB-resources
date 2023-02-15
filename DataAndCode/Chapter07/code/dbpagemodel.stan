data{
   int n;
   real x[n];
   real y[n];
   real mubeta0;
   real taubeta0;
   real mubeta1;
   real taubeta1;
   real a;
   real b;
}
transformed data{
   real sigmabeta0;
   real sigmabeta1;
   
   sigmabeta0 = 1/sqrt(taubeta0);
   sigmabeta1 = 1/sqrt(taubeta1);
}

parameters{
   real beta0;
   real beta1;
   real tau;
}
transformed parameters{
   real sigma;
   real mu[n];
   
   for (i in 1:n){
      mu[i] = beta0 + (beta1 * x[i]);
   }
   sigma = 1/sqrt(tau);
}

model{
  for(i in 1:n){
    y[i] ~ normal(mu[i], sigma);
  }
  beta0 ~ normal(mubeta0, sigmabeta0);
  beta1 ~ normal(mubeta1, sigmabeta1);
  tau ~ gamma(a,b);
}
