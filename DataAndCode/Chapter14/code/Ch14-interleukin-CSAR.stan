//  This program in Stan matches the code in JAGS and BUGS.
//  Running the program in Stan leads to some divergent transitions
//  of the sampler, however.
//  The specification of the correlations is problematic with these data.
//  There are more efficient or better ways to specify the program in Stan.
//  We have kept this version, however, for consistency with JAGS and BUGS.
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
    matrix[N,6] r;
    real med;
    real ten[2];
    real ninety[2];
    real rho1;
    real rho2;
    real sig2;
    real sigg2;
    real sigw2;

    for (i in 1:N){
       for (j in 1:6){
          r[i,j] = rhow; // If unequal time, multiply rhow by increments.
       }
    }
    sig2 = square(sig);
    sigg2 = square(sigg);
    sigw2 = square(sigw);
    rho1 = sigg2/(sig2 + sigg2);
    rho2 = (sigg2+(sigw2*rhow))/(sig2 + sigg2 + sigw2);
    ten[1] = beta1  - 1.28*sig;
    ninety[1] = beta1 + 1.28*sig;
    med = exp(beta1);
    ten[2] = exp(ten[1]);
    ninety[2] = exp(ninety[1]);
}
model {
   matrix[N,6] res;
   matrix[N,6] mu;
   
   beta1 ~ normal(0, 2);
   for (i in 1:N) {
      w[i,1] ~ normal(0, sigw);
      for (j in 2:6) {
         w[i,j] ~ normal((w[i,(j-1)]*r[i,(j-1)]), sigw*(1-square(r[i,(j-1)])));
      }
      gam[i] ~ normal(0,sigg);
      for (j in 1:6) {
         mu[i,j] = beta1 + gam[i] + w[i,j];
         logy[i,j] ~ normal(mu[i,j], sig);
	 res[i,j] = logy[i,j] - mu[i,j];
      }
   }    
}
