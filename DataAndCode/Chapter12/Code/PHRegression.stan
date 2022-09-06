data {
    int n;
    int k;
    real w;
    vector<lower=0>[n] y;
    vector<lower=0>[k] d;
    vector<lower=0>[n] AGgroup;
    int delta[n];
}
transformed data {
    vector[k+1] a;
    matrix[n,k] H;
    int N[n,k];
    matrix[n,k] min_mat;
    a[1] = 0;
    for (i in 2:k) {
      a[i] = (d[i-1] + d[i])/2;
    }
    a[k+1] = d[k] + (d[k] - d[k-1])/2;
    for (i in 1:n) {
      for (j in 1:k) {
          if (y[i] < a[j+1])
              min_mat[i,j] = y[i];
          else 
             min_mat[i,j] = a[j+1]; 
         //H[i,j] = HH[i,j]*(min(y[i],a[j+1])-a[j]);
         H[i,j] = (y[i] >= a[j])*(min_mat[i,j]-a[j]);
         N[i,j] = (y[i] >= a[j])*(a[j+1] > y[i])*delta[i];
      }
   }
}
parameters {
    real beta;
    vector<lower=0>[k] lam;
}
transformed parameters {
    vector[n] mu;
    real<lower=0> theta[n,k];
    for (i in 1:n) {
      mu[i] = beta * (AGgroup[i]);
      for (j in 1:k) {
         theta[i,j] = exp(mu[i])*H[i,j]*lam[j];
      }
    }
}
model {
    for (i in 1:n) {
      for (j in 1:k) {
          if(theta[i,j] > 0)
            N[i,j] ~ poisson(theta[i,j]);
          else
            target += 0;
      }
    }
    lam[1] ~ gamma(100,10000);
    for (j in 2:k) { 
       lam[j] ~ gamma(0.05*w*(a[j+1]-a[j]), w*(a[j+1]-a[j]));
    }
    beta ~ normal(0, 1/sqrt(0.000001));
}
generated quantities {
    real HR;
    HR = exp(beta);
}
