
Example 8.1
file to be named 

model{
for (i in 1:2) {LE[i] ~ dbin(p[i],n[i])
logit(p[i]) <- beta[1] + beta[2]*High[i] }
beta[1]  ~ dflat()
beta[2]  ~ dflat()
theta[2]<- exp(beta[1])/(1+exp(beta[1]))
theta[1]<- exp(beta[1]+beta[2])/(1+exp(beta[1]+beta[2]))
# Candy
OR <- exp(beta[2])
RR <- theta[1]/theta[2]
RD <- theta[1] - theta[2]
}
# Data
list(n = c(477,734),High = c(1,0), LE= c(126,53))
#Initial values for theta's
list(beta=c(0,0))
______________________________________________________________

Example 8.2
file to be named

model { 
# LR model for the data with sAge
# LE[i] = indicator of lymphedema for ith woman
# Nexam[i] = number of lymph notes examined
# Met[i] = indicator of metastasis in lymph node
# Age[i] = age of ith woman
  for (i in 1:1210)  {
    LE[i] ~ dbern(p[i])
    logit(p[i]) <- beta[1] +
    beta[2]*step(Nexam[i] -6)  # This dichotomizes Nexam
              + beta[3]*Met[i] +
    beta[4]*(Age[i] -mean(Age[]))/sd(Age[])   
  }

# Prior
  for (i in 1:4) {beta[i] ~ dflat()  } # Flat prior on beta

# 
  OR[1] <- exp(beta[2])  
  OR[2] <- exp(beta[3])
  OR[3] <- exp(beta[4])

# Probabilities of LE corresponding to ave Age, 72.7
  PrLE[1] <- exp(beta[1]+beta[2] + beta[3])/
              (1+ exp(beta[1] +beta[2]+ beta[3]))# (High,Met)
  PrLE[2] <- exp(beta[1] + beta[3])/
              (1+ exp(beta[1] + beta[3])) #M-L
  RR[1] <- PrLE[1]/PrLE[2] #Effect of H/L for Met
  PrLE[3] <- exp(beta[1]+beta[2])/
              (1+ exp(beta[1] +beta[2])) #(No Met,H)
  PrLE[4] <- exp(beta[1])/(1+ exp(beta[1]))# (No Met,Low)
  RR[2] <- PrLE[3]/PrLE[4] # Effect of H/L for No Met

# Probabilities of LE Corresponding to Age 78.3
  PrLE[5] <- exp(beta[1]+beta[2] + beta[3] + beta[4])/
              (1+ exp(beta[1] +beta[2] + beta[3] + beta[4]))#(Met,High)
  PrLE[6] <- exp(beta[1]+beta[3] + beta[4])/
              (1+ exp(beta[1] +beta[3]+ beta[4])) # (Met,Low)
  RR[3] <- PrLE[5]/PrLE[6] # Effect of H/L for Met

  PrLE[7] <- exp(beta[1]+beta[2] + beta[4])/
               (1+ exp(beta[1] +beta[2]+ beta[4])) # (No Met,High)
  PrLE[8] <- exp(beta[1] + beta[4])/
              (1+ exp(beta[1] + beta[4])) # (No Met,Low)
  RR[4] <- PrLE[7]/PrLE[8] # Effect of H/L for No-Met

  prob[1] <- step(beta[2])  
  prob[2] <- step(beta[3] )
  prob[3] <- step(beta[4])                  
}

list(beta = c(0,0,0,0))
LE[]Nexam[] Met[] Age[]
1	  24	0	67
0	   1	0	70
0	   1	1	75
0	   2	0	67
0	  11	0	73
0	   0	0	67
...
END

________________________________________________________

Example 8.5
file to be named 

model{
  for (i in 1:4) { 
    y[i] ~ dbin(p[i],n[i])
    logit(p[i]) <- beta[1] + beta[2]*sDose[i] 
  }
  for (i in 1:2) { 
    beta[i] ~ dunif(-10,10) 
  }
  ED50 <- (logit(0.5) - beta[1])/beta[2]

#  The following code is used to modify the model
# to use Uniform priors on the theta's
# Note that we assume theta[2] < theta[1]
# theta[1] <- expit( beta[1] ) 
# expit(a) = e^a/(1+e^a)
# theta[2] <- expit( beta[1] + beta[2]*2)
# beta[1]  <- logit(theta[1])
# beta[2] <- (logit(theta[2]) - beta[1])/2
# theta[1] ~ dunif(0,1)
# theta[2] ~ dunif(0,theta[1])
}
list(y = c(8,5,3,0), n = c(8,8,8,8), sDose = c(6,2,0,-1))
list(beta = c(0,0))
list(theta = c(0.5, 0.25)) # theta[2] < theta[1]
_____________________________________________________________

Example 8.6
file to be named 

model{ 
  for (i in 1:2) {
    theta[i] ~ dbeta(a[i],b[i])
  }
  beta0 <- logit(theta[1])
  beta1 <- logit(theta[2]) - beta0
}
list(a = c(2.044,2.044), b = c(13.006,13.006))

and these lines within outer braces {\tt model\{\}} can be included in the code to analyze the data (which will, of course, need to include the data sampling model, the data, and the initial values). We could also use the {\tt R} code

theta1 <- rbeta(5000,2.004,13.006)
theta2 <- rbeta(5000,2.004,13.006)
beta0 <- log(theta1/(1-theta1))
beta1 <- log(theta2/(1-theta2)) - beta0
______________________________________________________________

Example 8.7

\begin{verbatim}
# R code for generating Xtinv
Xtilde <- c(1,25,7.84,60,0,0,
            1,25,3.34,10,0,0,
            1,41,3.34,60,1,60,
            1,41,7.84,10,1,10,
            1,33,5.74,35,0,0,
            1,33,5.74,35,1,35)
Xtilde <- matrix(Xtilde, nrow=6, ncol=6, byrow=TRUE)
Xtinv <- solve(Xtilde)
\end{verbatim}

model{
 for(i in 1:n) {
    death[i] ~ dbern(theta[i])
    logit(theta[i]) <- beta[1] + beta[2]*ISS[i]
        + beta[3]*RTS[i] + beta[4]*AGE[i]
        + beta[5]*TI[i] + beta[6]*AGE[i]*TI[i]    
  }
  for(i in 1:6) {
    tildetheta[i] ~ dbeta(a[i],b[i])
    v[i] <- log(tildetheta[i]/(1-tildetheta[i]))
    beta[i] <- inprod(Xtinv[i,1:6], v[1:6])
  }
}  
list(tildetheta=c(0.5,0.5,0.5,0.5,0.5,0.5))
list(n=300, a=c(1.1,3,5.9,1.3,1.1,1.5), b=c(8.5,11,1.7,12,4.9,5.5))
Xtinv[,1] Xtinv[,2] Xtinv[,3] Xtinv[,4] Xtinv[,5] Xtinv[,6]
  4.647917  6.047917  3.285417  3.285417  -9.695833  -6.570833
 -0.031250 -0.031250  0.031250  0.031250   0.062500  -0.062500
 -1.666667 -1.666667 -1.666667 -1.666667   3.333333   3.333333
  0.170000  0.130000  0.150000  0.150000  -0.300000  -0.300000
 11.200000  9.800000  9.800000 11.200000 -22.000000 -20.000000
 -0.320000 -0.280000 -0.280000 -0.320000   0.600000   0.600000
END
ID[ ]   death[ ]   ISS[ ]   TI[ ]    RTS[ ]    AGE[ ]
2979      0          1       1       7.8408      25
1167      0          9       0       7.8408      63
116       0         29       0       2.9304      32
....
END

