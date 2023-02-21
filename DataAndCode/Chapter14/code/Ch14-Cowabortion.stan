data {
    int N;
    vector[N] GR;
    vector[N] DO;
    int<lower=0> y[N];
    int<lower=0> herd[N];
}
transformed data{
    vector[N] sGR;
    vector[N] sDO;
    sGR = (GR - 3.485024)/1.517102;
    sDO = (DO - 115.842)/59.72637;
    //sGR = (GR - mean(GR))/sd(GR);
    //sDO = (DO - mean(DO))/sd(DO);
}
parameters {
    vector[4] beta;
    real<lower=0,upper=2> sigmag;
    vector[9] gamma;
}
transformed parameters {
    vector<lower=0,upper=1>[N] theta;
    for(k in 1:N) {
      theta[k] =  inv_logit(beta[2]*sDO[k] + beta[3]*sGR[k] + beta[4]*sDO[k]*sGR[k] + gamma[herd[k]]);
    }
}
// This is the cow abortion data from chapter 14 of the book.
model {
   for (k in 1:9) { 
      gamma[k] ~ normal(beta[1], sigmag);
   }
   sigmag ~ uniform(0,2);
   for (k in 1:N) {  
      y[k] ~ bernoulli(theta[k]);
      
   }
   beta[1] ~ normal(-2,1);
   for (i in 2:4) { 
       beta[i] ~ normal(0,1);
   }
}
generated quantities {
    vector[4] prob;
    for (i in 1:4) { 
      prob[i] = beta[i] > 0;  
   }
}

