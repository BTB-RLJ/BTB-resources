//  Stan program for Model 4 in Table 10.2.  Obtains CPOinv, which is
//  used to obtain LPML in post processing the output. 
data{
   int n;
   int p;
   real a;
   real b;
   int Smoke[n];
   real FEV[n];
   real Age[n];
   vector[p] beta0;
   matrix[p,p] C0;
}

parameters{
   vector[p] beta;
   real tau;
}

transformed parameters{
   real sigma;
   real RM;
   real MD;
   real mu[n];
   real meanFEVs[10];
   real meanFEVns[10];
   
   sigma = 1/sqrt(tau);

   for (i in 1:n){
      mu[i] = beta[1] + beta[2]*Age[i] + beta[3]*Smoke[i]  + beta[4]*Age[i]*Smoke[i];
   }
   
// Estimate mean FEV for smokers and nonsmokers who are 10, ..., 19 years old   
   for (i in 1:10){
      meanFEVs[i] =  beta[1] + (beta[2]+beta[4])*(i+9) + beta[3];
      meanFEVns[i] = beta[1] + beta[2]*(i+9);
   }

// Estimate relative mean comparing 18 year old smoker to 18 year old nonsmoker
   RM = meanFEVns[9]/meanFEVs[9];
   
// Estimate mean difference of 18 year old smoker minus 13 year old nonsmoker
   MD = meanFEVs[9]-meanFEVns[4]; 
}

model{
   for (i in 1:n){
      FEV[i] ~ normal(mu[i], sigma);
   }
   beta ~ multi_normal(beta0, C0);
   tau ~ gamma(0.001, 0.001);
}

generated quantities{
   real mu20s;
   real mu20ns;
   real FEV20s;
   real FEV20ns;

// Predict the FEV for a 20 year old smoker and nonsmoker
   mu20s =  beta[1] + (beta[2]+beta[4])*20 + beta[3];
   mu20ns = beta[1] + beta[2]*20;

   FEV20s = normal_rng(mu20s, sigma);
   FEV20ns = normal_rng(mu20ns, sigma);
}
