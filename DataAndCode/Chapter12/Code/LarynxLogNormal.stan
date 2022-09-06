data {
    int<lower=0> Nobs;
    int<lower=0> Ncens;
    vector<lower=0>[Nobs] tobs;
    vector<lower=0>[Ncens] tcens;
    vector[90] age;
    vector[90] yr;
    int stageobs[Nobs];
    vector[Nobs] ageobs;
    vector[Nobs] yrobs;
    int stagecens[Ncens];
    vector[Ncens] agecens;
    vector[Ncens] yrcens;
}
transformed data{
    vector[Nobs] sAgeobs;
    vector[Nobs] sYrobs;
    vector[Ncens] sAgecens;
    vector[Ncens] sYrcens;
    sAgeobs = (ageobs-mean(age))/sd(age);
    sYrobs = (yrobs-mean(yr))/sd(yr);
    sAgecens = (agecens-mean(age))/sd(age);
    sYrcens = (yrcens-mean(yr))/sd(yr);
}
parameters {
    vector[6] beta;
    real<lower=0> tau;
}
transformed parameters {
    real sigma;
    vector[Nobs] muobs;
    vector[Ncens] mucens;
    sigma = sqrt(1/tau);
    for(i in 1:Nobs) {
         muobs[i] = beta[1] + beta[2] * (stageobs[i] == 2) + beta[3]*(stageobs[i] == 3) + beta[4]*(stageobs[i] == 4) + beta[5]*sAgeobs[i] + beta[6]*sYrobs[i];
    }
    for(i in 1:Ncens) {
         mucens[i] = beta[1] + beta[2] * (stagecens[i] == 2) + beta[3]*(stagecens[i] == 3) + beta[4]*(stagecens[i] == 4) + beta[5]*sAgecens[i] + beta[6]*sYrcens[i];
    }
}
model {
    for(i in 1:Nobs) {
       tobs[i] ~ lognormal(muobs[i], sigma); 
    }
    for(i in 1:Ncens) {
      target += lognormal_lccdf(tcens[i]| mucens[i], sigma);  
    }
    for(i in 1:6) {
       beta[i] ~ normal(0, 1 / sqrt(0.000001)); 
    }
    tau ~ gamma(0.001,0.001);
}
