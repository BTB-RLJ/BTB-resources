data{
   int n[2];
   int y[2];
   real a[2];
   real b[2];
}

parameters{
   real theta[2];
}
transformed parameters{
   real odds[2];
   real RD;
   real RR;
   real OR;
   real test;

   for (i in 1:2){
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
      theta[i] ~ beta(a[i], b[i]);
   }
}

   
