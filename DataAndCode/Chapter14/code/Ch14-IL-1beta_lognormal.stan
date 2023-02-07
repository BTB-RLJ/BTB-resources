data {
    int N;
    int nID;
    vector[N] y; // Assumes data NOT on log scale
    int ID[N];
}
parameters {
    vector[nID] gam;
    real<lower=0, upper=2> sig;
    real beta1; 
    real<lower=0, upper=2> sigg;
}
transformed parameters {
    real med;
    real rho;
    real taug;
    real tau;
    real ten[2];
    real ninety[2];
    
    med = exp(beta1);
    taug = 1/sigg ^ 2;
    tau = 1 / sig ^ 2;
    rho = (1/taug)/(1/tau + 1/taug);
    ten[1] = beta1  - 1.28*sig;
    ninety[1] = beta1 + 1.28*sig;
    ten[2] = exp(ten[1]);
    ninety[2] = exp(ninety[1]);
}
model {
   for (i in 1:N) {
//      logy[i] ~ normal(gam[ID[i]], sig);
      y[i] ~ lognormal(gam[ID[i]], sig);      
   }    
   for (i in 1:nID) {
      gam[i] ~ normal(beta1, sigg) ; 
   }
   beta1 ~ normal(0, 2);
   sig ~ uniform(0, 2);
   sigg ~ uniform(0, 2);
}
generated quantities {
   real gamf;
   real logyf;
   real yf;
   real muf;
   
   gamf = normal_rng(beta1, sigg);
//   logyf = normal_rng(gamf, sig);
//   yf = exp(logyf);
   yf = lognormal_rng(gamf, sig);
   muf = exp(gamf);
}
