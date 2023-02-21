data {
  int n;        // Number of observations
  vector[n] y;  // Dependent variable
  real alpha0;
  real lambda0;
  real mu0;
  real tau0;   // Prior precision of mu for consistency with WinBUGS
}

parameters {
  real mu;
  real<lower=0> precision;
}

transformed parameters{
// Stan uses the standard deviaion, not the variance or the precision.
	real sigma;
	real cv;
	
	sigma = 1/sqrt(precision);
	cv = sigma/mu;
}

model {
  mu ~ normal(mu0, sqrt(1/tau0));
  precision ~ gamma(alpha0, lambda0);
  for (i in 1:n) {
    y[i] ~ normal(mu, sigma);
  }
}
