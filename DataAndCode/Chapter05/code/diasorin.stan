
data{
   int nL;
   int nN;
   real bL;
   real bN;
   real c; // c, d, e, f are parameters in gamma priors for precisions
   real d;
   real e;
   real f;
   real LowTurnover[nL];
   real NormalTurnover[nN];
}

transformed data{
   real logLowTurnover[nL];
   real logNormalTurnover[nL];

   for(i in 1:nL){
      logLowTurnover[i] = log(LowTurnover[i]);
   }
  for(i in 1:nN){
      logNormalTurnover[i] = log(NormalTurnover[i]);
  }
}

parameters{
   real muL;
   real muN;
   real tauL;
   real tauN;
}

transformed parameters{
   real Delta;
   real medL;
   real medN;
   real relmed;
   real sigmaL;
   real sigmaN;
   real rsd;
   real test[2];

   Delta = muN - muL;
   medL = exp(muL);          
   medN = exp(muN);
   relmed = medN /medL;        

   sigmaN = sqrt(1/tauN);     
   sigmaL = sqrt(1/tauL);
   rsd = sigmaL/sigmaN;
   test[1] = (rsd > 1.0) ? 1:0; // test[1] = step(rsd -1)
   test[2] = (relmed>1.0) ? 1:0; // test[2] = step(relmed -1)
}

model {
   for (i in 1:nL){
      logLowTurnover[i] ~ normal(muL, sigmaL);
   }
   for (i in 1:nN){
      logNormalTurnover[i] ~ normal(muN, sigmaN);
   }
   muL ~ normal(4.87, bL);  
   muN ~ normal(5.39, bN);
   tauL ~ gamma(c, d);   
   tauN ~ gamma(e, f);
}

generated quantities{
   real pdLlog;
   real pdNlog;
   real pdL;
   real pdN;
   real AUC;
   
 //  Sample predictive distributions Low and Normal
   pdLlog = normal_rng(muL, sigmaL); // on log-measurement axis     
   pdNlog = normal_rng(muN, sigmaN); // on log-measurement axis
   pdL = exp(pdLlog); //lognormal_rng(muL, sigmaL); // on original measurement axis     
   pdN = exp(pdNlog); //lognormal_rng(muN, sigmaN); // on original measurement axis
   AUC = (pdN > pdL) ? 1:0; // returns 1 if pdN > pdL, 0 otherwise.
}
