data {
    int<lower=0> n[2];
    real<lower=0> a1;
    real<lower=0> a2;
    real<lower=0> b1;
    real<lower=0> b2;
    vector<lower=0>[n[1]] t1;
    vector<lower=0>[n[2]] t2;
}
parameters {
    real<lower=0> theta1;
    real<lower=0> theta2;
}
model {
    for(i in 1:n[1]) {
        t1[i] ~ exponential(theta1);
    }
    for(i in 1:n[2]) {
        t2[i] ~ exponential(theta2);
    }
    theta1 ~ gamma(a1,b1);
    theta2 ~ gamma(a2,b2);
}
generated quantities {
   real median1;
   real median2;
   real relmedian;
   vector[2] S;
   real Sdiff;
   median1 = log(2)/theta1;
   median2 = log(2)/theta2;
   relmedian = theta2/theta1;
   S[1] = exp(-24*theta1);
   S[2] = exp(-24*theta2);
   Sdiff = S[1]-S[2];
}