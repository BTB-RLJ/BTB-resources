data{
	int n[2];
	vector[n[1]+n[2]] y;
	vector[n[1]+n[2]] age;
}
transformed data{
	int N;
	real meanAge;
	real sdAge;
	vector[n[1]+n[2]] sAge;
	N = n[1]+n[2];
	meanAge = mean(age);
	sdAge = sd(age);
	for (i in 1:N) {
		sAge[i] = (age[i]-meanAge)/sdAge;
	}
}
parameters{
	real b0[2];
	real b1[2];
	vector<lower=0>[2] tau;
}
transformed parameters{
	vector<lower=0>[2] sigma;
	real prob0;
	real prob1;
	real mu0[n[1]];
	real mu1[n[2]];
	real b;
    	for(i in 1:2) {
       	      sigma[i] = 1 / sqrt(tau[i]); 
    	}
	b = sigma[1]/sigma[2];
    	for (i in 1:n[1]){
    	    mu0[i] = b0[1] + b0[2]*sAge[i];
    	}
   	for (i in (n[1] +1):(n[1] + n[2])){
//          mu1[(i-n[1])] = b1[1];           // Model 1 removes the slope     
      	    mu1[(i-n[1])] = b1[1] + b1[2]*sAge[i]; // Model 2 includes the slope
    	}
	prob0 = (b0[2]>0);
	prob1 = (b1[2]>0);  
}

model{
	for (i in 1:n[1]){
		y[i] ~ normal(mu0[i], sigma[1]);
	}  
	for (i in (n[1] +1):(n[1] + n[2])){
		y[i] ~ normal(mu1[(i-n[1])], sigma[2]) ;
	}
// SPECIFY PRIORS
	b0[1] ~ normal(0, 1000);
	b0[2] ~ normal(0, 1000);
	b1[1] ~ normal(0, 1000);
	b1[2] ~ normal(0, 1000);
	tau[1] ~ gamma(.001,.001);
	tau[2] ~ gamma(.001,.001);
}
generated quantities{
	real x[2];
	real a[2];
	real c[2];
	real AUC[2];
	real diffAUC;
	real probAUC;
// Obtain AUC for two ages
	x[1] = 70;
	x[2] = 40;
	for (i in 1:2) {
// For model 1 that does not include a slope for group 2
//		a[i] =( (b1[1] - b0[1]) +  ( (0.0-b0[2])*(x[i] -meanAge )/sdAge) )/sigma[2];
// For model 2 that includes a slope for group 2
		a[i] =( (b1[1] - b0[1]) +  ( (b1[2]-b0[2])*(x[i] -meanAge )/sdAge) )/sigma[2]; 
		c[i] =  a[i]/sqrt(1 + (b*b));
		AUC[i] = normal_cdf(c[i], 0, 1);
	}     
	diffAUC = AUC[2] - AUC[1];  
	probAUC = (diffAUC > 0);
}
