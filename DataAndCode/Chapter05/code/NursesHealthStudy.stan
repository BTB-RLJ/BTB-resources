data{
   int y[2];
   real M[2];
}

parameters{
   real theta[2];
}
transformed parameters{
   real theta1;
   real theta2;
   real RateRatio;
   real test;
   
   theta1 = theta[1]*M[1];
   theta2 = theta[2]*M[2];
   RateRatio = theta[1]/theta[2];
   test =  (RateRatio>1) ? 1:0;
}

model {
   y[1] ~ poisson(theta1);
   y[2] ~ poisson(theta2);
   theta[1] ~ gamma(0.001,0.001);
   theta[2] ~ gamma(0.001,0.001);
 }
