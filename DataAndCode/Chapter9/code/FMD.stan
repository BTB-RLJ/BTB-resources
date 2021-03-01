data {
    vector[66] Size;
    int<lower=0> FMD[66];
    vector[66] East;
}
transformed data {
    vector[66] sSize;
    sSize = (Size - mean(Size))/sd(Size);
}
parameters {
    vector[3] b;
}
transformed parameters {
    vector[66] lambda;
    for(i in 1:66) {
        lambda[i] = exp(b[1] + b[2]*East[i] + b[3]*sSize[i]);
    }
}
model {
    for(i in 1:66) {
        FMD[i] ~ poisson(lambda[i]);
    }
    for (j in 1:3) { 
      b[j] ~ uniform(-2,2);
    }
}
generated quantities {
    vector[66] pp;
    vector[66] q;
    vector[66] o;
    real Pearson;
    real Ezeros;
    real Obszeros;
    for(i in 1:66) {
      pp[i] = ((FMD[i] - lambda[i])/sqrt(lambda[i]))^2;  
    }
    
    for(i in 1:66) {
        q[i] = exp(-lambda[i]);
    }
    
    for(i in 1:66) {
       o[i] = FMD[i] == 0;
    }
    Pearson = sum(pp); 
    Ezeros = sum(q);
    Obszeros = sum(o); 
}