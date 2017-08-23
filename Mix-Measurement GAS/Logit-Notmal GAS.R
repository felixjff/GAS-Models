#####Purpose: Simulation and calibration of a mixed-measurement observation driven GAS model to determine the bias of MLE estimates
###Author: Felix Farias Fueyo
###Details: 
# 1. The mix-measurement model consists of a normal and a logit component.
# 2. The logit and normal components are drived by a shared factor

###              ###
###  Simulation  ###
###              ###

simulations <- 1000
tobs <- 100
cobs <- 100

#Create storage space for the scores, factor, normal variable, success indicator, true success rate and generated success rate
#Common
score <- matrix(data = NA, nrow = tobs, ncol = simulations)
f <- matrix(data = NA, nrow = tobs+1, ncol = simulations)
#Logit
score_l <- matrix(data = NA, nrow = tobs, ncol = simulations)
p <- matrix(data = NA, nrow = tobs, ncol = simulations)
gen_success_rate <- matrix(data = NA, nrow = tobs, ncol = simulations)
#Normal
score_n <- matrix(data = NA, nrow = tobs, ncol = simulations)
y_n <- matrix(data = NA, nrow = tobs ,  ncol = simulations)

# Data Generating Process True Parameter Values
true_param <- c(f1 = 0 , Zm = 0.5, Zc = 0.5, w = 0, A = 0.5, B = 0.5, Sigma_squared = 1)
f1 = 0 ; Zm = 0.5; w = 0; A = 0.5; B = 0.5; Sigma_squared = 1; Zc = 0.5

#simulate
for(i in 1:simulations){
  #initialize factor
  f[1,i] <- f1
  
  #generate path
  for(t in 1:tobs){
    #Dynamic probability for logit component
    p[t,i] <- 1/(1+ exp(-Zc*f[t,i]))
    
    #Instead of storing all observations, only store the generated sucess rate.
    y_l <- rbinom(n = cobs, size = 1, prob = p[t,i])
    gen_success_rate[t,i] <- sum(y_l)/cobs
    
    #Draw from normal distribution with dynamic mean specification. 
    y_n[t,i] <- rnorm(1, mean = Zm*f[t,i], sd = Sigma_squared)
    
    #Compute scores
    score_l[t,i] <- sum(y_l)*Zc - cobs*p[t,i]*Zc
    score_n[t,i] <- (1/Sigma_squared)*Zm*(y_n[t,i] - Zm*f[t,i])
    
    #Compute score of mixed-measurement log-likelihood function
    score[t,i] <- score_n[t,i] + score_l[t,i]
    
    
    f[t+1,i] <- w + A*score[t,i] + B*f[t,i]
  }
  print(i)
}

#Note how the simulation results support the theoretical results, namely
## 1. The score series for both distributions are martingales with mean zero. 
## 2. The unconditional mean of the factor is E[f] = w/(1-B)



###              ###
###  Estimation  ###
###              ###

#Define the loglikelihood function for the systematic-logit GAS model 
loglikelihood <- function(par, path_l, path_n, tobs, cobs){
  f1 <- par[grepl("f1",names(par))]
  Zc <- par[grepl("Zc",names(par))]
  Zm <- par[grepl("Zm",names(par))]
  w <- par[grepl("w",names(par))]
  A <- par[grepl("A",names(par))]
  B <- par[grepl("B",names(par))]
  ssq <- par[grepl("Sigma_squared",names(par))]

  #Common
  score_ <- matrix(data = NA, nrow = tobs, ncol = 1)
  f_ <- matrix(data = NA, nrow = tobs+1, ncol = 1)
  loglike <- matrix(data = NA, nrow = tobs, ncol = 1)
  #Logit
  p_ <- matrix(data = NA, nrow = tobs, ncol = 1)
  score_l_ <- matrix(data = NA, nrow = tobs, ncol = 1)
  #Normal
  y_n_ <- matrix(data = NA, nrow = tobs ,  ncol = simulations)
  score_n_ <- matrix(data = NA, nrow = tobs, ncol = 1)
  
  #Initialize GAS component
  f_[1] <- f1
  
  #compute likelihood and other elements at every t
  for(i in 1:tobs){
    #Dynamic probability for logit component
    p_[i] <- 1/(1 + exp(-Zc*f_[i]))
    
    #Score
    score_l_[i] <- cobs*path_l[i]*Zc - cobs*p_[i]*Zc
    score_n_[i] <- (1/ssq)*Zm*(path_n[i] - Zm*f_[i])
    
    score_[i] <-  score_l_[i] + score_n_[i]
    
    #Log-likelihood
    loglike_l <- cobs*path_l[i]*Zc*f_[i] - cobs*log(1 + exp(Zc*f_[i]))
    loglike_n <- -0.5*log(2*pi*ssq) - 0.5*(1/ssq)*(path_n[i] - Zm*f_[i])^2
    
    loglike[i] <- loglike_n + loglike_l
    
    f_[i+1] <- w + A*score_[i] + B*f_[i]
  }
  
  #compute log-likelihood
  loglike <- sum(loglike)
  return(-loglike)
}

#Obtain estimation results for every path
param_0 <- c(f1 = 0.5 , Zm = 2, Zc = 1, w = 1, A = 0.1, B = 0.99, Sigma_squared = 3)

param_ests <- matrix(data = NA, nrow = simulations, ncol = length(param_0))

for(i in 1:simulations){
  fit <- optim(par = param_0, fn = loglikelihood, method = "BFGS", path_l =  gen_success_rate[,i],
               path_n = y_n[,i], tobs = tobs, cobs = cobs)
  param_ests[i,] <- fit$par 
  print(i)
}


#Compute bias for each parameter
param_bias <- colMeans(param_ests) - true_param
