data {
  int<lower=0> N; 
  real x[N]; 
  real Y[N]; 
} 
parameters {
  real alpha; 
  real beta;  
  real<lower=.5,upper= 1> gamma; 
  real<lower=0> tau; 
} 
transformed parameters {
  real sigma; 
  real U3; 
  sigma = 1 / sqrt(tau); 
  U3 = logit(gamma);
} 
model {
  real m[N];
  for (i in 1:N) {
     m[i] = alpha - beta * pow(gamma, x[i]); 
  } 
  Y ~ normal(m, sigma); 
  alpha ~ normal(0.0, 1000); 
  beta ~ normal(0.0, 1000); 
  gamma ~ uniform(.5, 1); 
  tau ~ gamma(.0001, .0001); 
}

