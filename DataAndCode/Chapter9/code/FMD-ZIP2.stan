// Zero-inflated Poisson
data {
    int<lower=0> N;
    int<lower=0> FMD[N];
    vector[N] East;
    vector[N] sSize;
}
transformed data{
    int<lower=0, upper=1> z[N];
    for (i in 1:N) {
        z[i] = (FMD[i] == 0.0);
    }
}
parameters {
    vector<lower=-5, upper=5>[3] b;
    vector<lower=-5, upper=5>[3] g;
}
transformed parameters {
    vector[N] lambda;
    vector[N] pi;
    for(i in 1:N) {
        lambda[i] = exp(b[1] + (b[2]*East[i]) + (b[3]*sSize[i]));
	pi[i] = inv_logit(g[1] + (g[2]*East[i]) + (g[3]*sSize[i]));
    }
}
model {
  for (i in 1:N) {
    if (FMD[i] == 0)
      target += log_sum_exp(bernoulli_lpmf(1 | pi[i]),
                            bernoulli_lpmf(0 | pi[i])
                              + poisson_lpmf(FMD[i] | lambda[i]));
    else
      target += bernoulli_lpmf(0 | pi[i])
                  + poisson_lpmf(FMD[i] | lambda[i]);
  }
}
generated quantities {
    vector[N] pp;
    vector[N] E;
    vector[N] PR;
    vector[3] bProbPos;
    vector[3] gProbPos;
    real Pearson;
    for(i in 1:N) {
      E[i] = z[i]*pi[i] + (1-z[i])*lambda[i];
      PR[i] = (FMD[i] - E[i])/sqrt(E[i]);
      pp[i] = PR[i]^2;
    } 
    for(i in 1:3) {
        bProbPos[i] = b[i]>0;
	gProbPos[i] = g[i]>0;
    }
    Pearson = sum(pp); 
}
