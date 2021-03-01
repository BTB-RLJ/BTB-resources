data {
    int N;
    vector[N] T0;
    real meanT0;
    vector[N] Y;
    real<lower=0> nu;
}
parameters {
    real b0;
    real b1;
    vector<lower=0>[N] lambda;
    real<lower=0> sigma2;
}
transformed parameters {
    vector<lower=0>[N] tau;
    vector[N] mu;
    for(i in 1:N) {
        tau[i] = 1 / (sigma2 / lambda[i]);
        mu[i] = b0 + b1*(T0[i] - meanT0);
    }
}
model {
    for (i in 1:N) {
		Y[i] ~ normal(mu[i], 1 / sqrt(tau[i]));
		lambda[i] ~ gamma(nu/2, 0.5);
	}
	b0 ~ normal(0, sqrt(10));
	b1 ~ normal(0, sqrt(10));
	sigma2 ~ inv_gamma(0.001, 0.001);
}
generated quantities {
    real StDev;
    StDev = sqrt(sigma2);
}