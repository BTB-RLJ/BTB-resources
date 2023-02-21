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
    matrix[N,6] gamm;
    real med;
    real rho1;
    real rho2;
    real sig2;
    real sigg2;
    real sigw2;
    real sigww2;
    real ten[2];
    real ninety[2];
    matrix[N,5] ww;
    
    for (i in 1:N){
       for (j in 1:6){
          gamm[i,j] = gam[i] + w[i,j];
       }
    }
    sig2 = sig^2;
    sigg2 = sigg^2;
    sigw2 = sigw^2;
    sigww2 = sigw2*(1-(rhow^2));
    rho1 = sigg2/(sig2 + sigg2);
    rho2 = (sigg2+(sigw2*rhow))/(sig2 + sigg2 + sigw2);
    ten[1] = beta1  - 1.28*sig;
    ninety[1] = beta1 + 1.28*sig;
    med = exp(beta1);
    ten[2] = exp(ten[1]);
    ninety[2] = exp(ninety[1]);
    for (i in 1:N) {
        for (j in 1:5) {
            ww[i,j] = rhow*w[i,j];
	}
    }
}
model {
   matrix[N,6] res;
//   matrix[N,5] ww;
   
   for (i in 1:N) {
      gam[i] ~ normal(beta1, sigg);
      logy[i,1] ~ normal(gamm[i,1], sig);
      res[i,1] = logy[i,1] - gamm[i,1];
      w[i,1] ~ normal(0, sigw);
      for (j in 2:6) {
         logy[i,j] ~ normal(gamm[i,j], sig);
	 res[i,j] = logy[i,j] - gamm[i,j];
//         ww[i,(j-1)] = rhow*w[i,(j-1)];
         w[i, j] ~ normal(ww[i,(j-1)], sqrt(sigww2));
      }
   }    
   beta1 ~ normal(0, 2);
}
