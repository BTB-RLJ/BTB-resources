data{
   int n;
   int y[4];
   real<lower=0> aSe;
   real<lower=0> bSe;
   real<lower=0> aSp;
   real<lower=0> bSp;
   real<lower=0> aPi;
   real<lower=0> bPi;   
}
parameters{
   real Se;
   real Sp;
   real pi;
}
transformed parameters{
   vector[4] p;
   p[1] = Se*pi;
   p[2] = (1-Sp)*(1-pi);
   p[3] = (1-Se)*pi;
   p[4] = Sp*(1-pi);
}
model{
   y ~ multinomial(p);
   Se ~ beta(aSe, bSe);
   Sp ~ beta(aSp, bSp);
   pi ~ beta(aPi, bPi);
}

