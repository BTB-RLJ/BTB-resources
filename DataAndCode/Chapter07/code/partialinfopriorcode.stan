data {
   int N; // Number of observations
   int p; // Number of covariates + 1 (for intercept)
   vector[N] Y;  // Dependent variable
   vector[2] m0; // Elicited expected values
}
parameters {
   vector[(p-2)] beta;
   vector[2] tildem;
}
transformed parameters {
   vector[2] beta12;
   beta12[1] <- tildem[1];
   beta12[2] <- tildem[2] - tildem[1];
}
model { 
   for (i in 1:2) { 
      tildem[i] ~ dnorm(m0[i],tau0[i]) 
   }
   for (j in 3:p) {
      beta[j] ~ dnorm(0,0.000001)
   }
//   ...  // Steps for modeling Y
}

// list(m0 = c(..),tau0= c(..), p = ) # Required "data" values
// list(tildem = c(...)) # Initial values for tildem
