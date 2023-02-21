data {
    int N;
    matrix[N, 6] logy;
}
parameters {
    real beta1; 
    real<lower=0, upper=1> rhow;
    real<lower=0, upper=2> sig;
    real<lower=0, upper=2> sigg;
    real<lower=0, upper=2> sigw;
    vector[N] gam;
    matrix[N,6] w;
}
transformed parameters {
    matrix[N,6] gamm;
    real med;
    real precw;
    real rho1;
    real rho2;
    real tau;
    real taug;
    real tauw;
    real ten[2];
    real ninety[2];
    
    for (i in 1:N){
       for (j in 1:6){
          gamm[i,j] = gam[i] + w[i,j];
       }
    }
    tau = (1/sig)^2;
    taug = (1/sigg)^2;
    tauw = (1/sigw)^2;
    precw = tauw/(1-(rhow^2));
    rho1 = (1/taug)/(1/tau + 1/taug);
    rho2 = ((1/taug)+((1/tauw)*rhow))/((1/tau) + (1/taug) + (1/tauw));
    ten[1] = beta1  - 1.28*sig;
    ninety[1] = beta1 + 1.28*sig;
    med = exp(beta1);
    ten[2] = exp(ten[1]);
    ninety[2] = exp(ninety[1]);
}
model {
   matrix[N,6] res;
   matrix[N,5] ww;
   
   for (i in 1:N) {
      gam[i] ~ normal(beta1, sigg);
      logy[i,1] ~ normal(gamm[i,1], sig);
      res[i,1] = logy[i,1] - gamm[i,1];
      w[i,1] ~ normal(0, sigw);
      for (j in 2:6) {
         logy[i,j] ~ normal(gamm[i,j], sig);
	 res[i,j] = logy[i,j] - gamm[i,j];
         ww[i,(j-1)] = rhow*w[i,(j-1)];
         w[i, j] ~ normal(ww[i,(j-1)], sqrt(1/precw));
      }
   }    
   beta1 ~ normal(0, 2);
}
