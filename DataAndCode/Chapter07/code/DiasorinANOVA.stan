data {
    int<lower=0> n[3];
    vector[n[1]] low;
    vector[n[2]] normal;
    vector[n[3]] high;
}
parameters {
    vector[3] mu;
    real<lower=0> tau;
}
transformed parameters {
    real P;
    real diff21;
    real diff31;
    real diff32;
    real prob21;
    real prob31;
    real prob32;
    real<lower=0>sigma;
    vector<lower=0>[3] med;
    
    for(i in 1:3) {
        med[i] = exp(mu[i]);
    }
    sigma = 1/sqrt(tau);
    diff21 = mu[2]-mu[1];
    diff31 = mu[3]-mu[1];
    diff32 = mu[3]-mu[2];
    prob21 = (diff21<0) ? 0:1; // ifelse(diff21 < 0, 0 ,1)
    prob31 = (diff31<0) ? 0:1; // ifelse(diff31 < 0, 0 ,1)
    prob32 = (diff32<0) ? 0:1; // ifelse(diff32 < 0, 0 ,1)
//  Test for ordering of medians
//    P = step(mu[3]-mu[2])*step(mu[2]-mu[1])
    P = ((mu[3]>mu[2]) && (mu[2]>mu[1])) ? 1:0;
}
model {
// Likelihood (lognormal model)
    for(i in 1:n[1]) {
        low[i] ~ lognormal(mu[1],sigma);
    }
    for(i in 1:n[2]) {
        normal[i] ~ lognormal(mu[2],sigma);
    }
    for(i in 1:n[3]) {
        high[i] ~ lognormal(mu[3],sigma);
    }
    // Prior
    mu[1]~normal(4.87, 1/sqrt(347.12));
    mu[2]~normal(5.39, 1/sqrt(357.42));
    mu[3]~normal(6.40, 1/sqrt(166.67));
    tau~gamma(0.001, 0.001);
}
generated quantities {
    real lowf_ls;
    real normalf_ls;
    real highf_ls;
    real lowf;
    real normalf;
    real highf;
    
    lowf_ls = normal_rng(mu[1], sigma);
    normalf_ls = normal_rng(mu[2], sigma);
    highf_ls = normal_rng(mu[3], sigma);
    lowf = exp(lowf_ls);
    normalf = exp(normalf_ls);
    highf = exp(highf_ls);
}
