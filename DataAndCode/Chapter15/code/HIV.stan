data {
    int<lower=0> y[3];
    real<lower=0> a;
    real<lower=0> b;
    real<lower=0> a1;
    real<lower=0> b1;
    real<lower=0> a2;
    real<lower=0> b2;
}
parameters {
    real<lower=0,upper=1> pi;
    real<lower=0,upper=1> se;
    real<lower=0,upper=1> sp;
}
transformed parameters {
    vector<lower=0,upper=1>[3] p;
    p[1] = pi*se;
    p[2] = (1-pi)*(1-sp);
    p[3] =  pi*(1-se) + (1-pi)*sp;
}
model {
    y ~ multinomial(p);
    pi ~ beta(a,b);
    se ~ beta(a1,b1);
    sp ~ beta(a2,b2);
}
generated quantities {
    real ppv;
    real pdgneg;
    ppv = p[1]/(1-p[3]);
    pdgneg = pi*(1-se)/p[3];
}