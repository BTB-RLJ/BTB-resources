data {
    int<lower=0> FMD[66];
//    int z[66];
    vector[66] East;
    vector[66] sSize;
}
transformed data{
    int<lower=0, upper=1> z[66];
    for (i in 1:66) {
        z[i] = (FMD[i] == 0.0);
    }
}
parameters {
    vector<lower=-2, upper=2>[3] b;
    vector<lower=-2, upper=2>[3] g;
}
transformed parameters {
    vector[66] lambda;
    vector[66] llambda;
    vector[66] pi;
    for(i in 1:66) {
        lambda[i] = exp(b[1] + (b[2]*East[i]) + (b[3]*sSize[i]));
	llambda[i] = ((1-z[i])*0.0) + (z[i]*lambda[i]);
	pi[i] = inv_logit(g[1] + (g[2]*East[i]) + (g[3]*sSize[i]));
    }
}
model {
    for(i in 1:66) {
	z[i] ~ bernoulli_logit(g[1] + (g[2]*East[i]) + (g[3]*sSize[i]));
        FMD[i] ~ poisson(llambda[i]);
    }
}
