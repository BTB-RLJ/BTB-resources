data {
    int nObs;
    int nID;
    int ID[nObs];
    vector[nObs] Age;
    vector[nObs] dist;
    vector[nObs] Male;
}
transformed data{
    vector[nObs] sAge;
    for (i in 1:nObs){
        sAge[i] = (Age[i] - 11)/2; // standardize Age;
    }
}
parameters {
    vector[3] beta;
    matrix[nID, 2] gam;
//    cholesky_factor_corr[2] Lcorr;
    vector<lower=0, upper = 10>[2] sigg;
    real<lower=-1,upper=1> rhog;
    real<lower=0, upper = 10> sig;
}
transformed parameters {
    vector[nObs] mu;
    matrix[2,2] Sigma;
//    vector[2] mu_g;
//    mu_g[1] = beta[1];
//    mu_g[2] = beta[2];
//    corr_matrix[2] R;  // correlation matrix
//    cov_matrix[2] Sigma;
//    R = multiply_lower_tri_self_transpose(Lcorr);
//    Sigma = quad_form_diag(R, sigg);
    Sigma[1,1] = sigg[1] ^ 2;
    Sigma[2,2] = sigg[2] ^ 2;
    Sigma[1,2] = rhog * sigg[1] * sigg[2];
    Sigma[2,1] = Sigma[1,2];
    for(i in 1:nObs) {
        mu[i] = gam[ID[i],1] + gam[ID[i],2]*sAge[i] + beta[3]*Male[i];
//        mu[i] = gam[ID[i],1] + beta[2]*sAge[i] + gam[ID[i],2]*sAge[i] + beta[3]*Male[i];
    }
}
model {
    for(i in 1:nObs) {
      dist[i] ~ normal(mu[i], sig);  
    }
    for(i in 1:nID) {
        gam[i,] ~ multi_normal(beta[1:2], Sigma);
    }
    for(i in 1:3) {
        beta[i] ~ normal(0, 10);
    }
    for(j in 1:2) {
        sigg[j] ~ uniform(0, 10);
    }
//    Lcorr ~ lkj_corr_cholesky(2.0); // prior for cholesky factor of correlation
    sig ~ uniform(0, 10);
    rhog ~ uniform(-1, 1);
}
