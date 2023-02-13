data{
   int x[2];
   int n[2];
   real mu;
   real prec;
   real a;
   real b;
}
transformed data{
   real sigma;

   sigma = 1/sqrt(prec);
}

parameters{
   real ttheta2; // Separate theta parameters, since prior only on ttheta[2];
   real delta;

}
transformed parameters{
   real ttheta[2];
   real OR;
   real O[2];

   ttheta[1] = exp(delta)*ttheta2/(1-ttheta2*(1-exp(delta)));
   ttheta[2] = ttheta2;
   O[1] = ttheta[1]/(1-ttheta[1]);
   O[2] = ttheta[2]/(1-ttheta[2]);
   OR = O[1]/O[2];
}

model{
   for (i in 1:2) {
      x[i] ~ binomial(n[i], ttheta[i]);
   }
//   x[1] ~ binomial(n[1], ttheta1);
//   x[2] ~ binomial(n[2], ttheta2);
   ttheta2 ~ beta(a,b);
   delta ~ normal(mu, sigma);
}
