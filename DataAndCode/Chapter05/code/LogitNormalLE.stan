data{
   int n[2];
   int y[2];
   real d;
}
transformed data{
   real sigma;

   sigma = 1/sqrt(d);
}

parameters{
   real beta[2];
}

transformed parameters{
   real theta[2];
   real odds[2];
   real RD;
   real RR;
   real OR;
   real test;

   for (i in 1:2){
      theta[i] = exp(beta[i])/(1+ exp(beta[i]));
      odds[i] = theta[i]/(1-theta[i]);
   }
   RD = theta[1] - theta[2]; // High exposure - low exposure
   RR = theta[1]/theta[2];
   OR = odds[1]/odds[2];

   test = (RD>0) ? 1:0;
}

model{
   for (i in 1:2){
      y[i] ~ binomial(n[i], theta[i]);
      beta[i] ~ normal(0, sigma);
   }
}

   
