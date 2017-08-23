#####Purpose: Simulation and calibration of a systematic-logit GAS model to determine the bias of MLE estimates
###Author: Felix Farias Fueyo


###              ###
###  Simulation  ###
###              ###

#Parameters of true data generating process (DGP)
true_param <- c(f1 = 0.5, w = 0.3, A = 0.01, B = 0.5, Zc = 0.5)
f1 = 0.5; w = 0.3; A = 0.01; B = 0.5; Zc = 0.5

#simulations, cross-section and time-series length
simulations <- 1000
cobs <- 100
tobs <- 100

#Define space for score, factors and realizations for the Binomial indicator variable
score <- matrix(data = NA, nrow = tobs, ncol = simulations)
f <- matrix(data = NA, nrow = tobs+1, ncol = simulations)
p <- matrix(data = NA, nrow = tobs, ncol = simulations)
gen_success_rate <- matrix(data = NA, nrow = tobs, ncol = simulations)

#simulate
for(i in 1:simulations){
  #initialize factor
  f[1,i] <- f1
  
  #generate path
  for(t in 1:tobs){
    p[t,i] <- 1/(1+ exp(-Zc*f[t,i]))
    
    y <- rbinom(n = cobs, size = 1, prob = p[t,i])
    
    #Instead of storing all observations, only store the generated sucess rate.
    gen_success_rate[t,i] <- sum(y)/cobs
    
    score[t,i] <- sum(y)*Zc - cobs*p[t,i]*Zc
    f[t+1,i] <- w + A*score[t,i] + B*f[t,i]
  }
  print(i)
}

#Note how the simulation results support the theoretical results, namely
## 1. The score series is a martingale difference with mean zero. 
## 2. The unconditional mean of the factor is E[f] = w/(1-B)


###              ###
###  Estimation  ###
###              ###

#Define the loglikelihood function for the systematic-logit GAS model 
loglikelihood <- function(par, path, tobs, cobs){
  f1 <- par[grepl("f1",names(par))]
  Zc <- par[grepl("Zc",names(par))]
  w <- par[grepl("w",names(par))]
  A <- par[grepl("A",names(par))]
  B <- par[grepl("B",names(par))]
  
  #create storage space for the score, GAS component and the log likelihood
  score_ <- matrix(data = NA, nrow = tobs, ncol = 1)
  f_ <- matrix(data = NA, nrow = tobs+1, ncol = 1)
  p_ <- matrix(data = NA, nrow = tobs, ncol = 1)
  loglike <- matrix(data = NA, nrow = tobs, ncol = 1)
  
  #Initialize GAS component
  f_[1] <- f1
  
  #compute likelihood and other elements at every t
  for(i in 1:tobs){
    p_[i] <- 1/(1 + exp(-Zc*f_[i]))
    
    score_[i] <- cobs*path[i]*Zc - cobs*p_[i]*Zc
    
    loglike[i] <- cobs*path[i]*Zc*f_[i] - cobs*log(1 + exp(Zc*f_[i]))
    
    f_[i+1] <- w + A*score_[i] + B*f_[i]
  }
  
  #compute log-likelihood
  loglike <- sum(loglike)
  return(-loglike)
}

#Obtain estimation results for every path
#param_0 <- c(f1 = 0, w = 0, A = 0.5, B = 0.99, Zc = 0.5)
param_0 <- c(f1 = 0.5, w = 0.3, A = 0.01, B = 0.5, Zc = 0.5)

param_ests <- matrix(data = NA, nrow = simulations, ncol = length(param_0))

for(i in 1:simulations){
  fit <- optim(par = param_0, fn = loglikelihood, method = "BFGS", path = gen_success_rate[,i], tobs = tobs, cobs = cobs)
  param_ests[i,] <- fit$par 
  print(i)
}

par <- c(f1 = 0, w = colMeans(param_ests, na.rm = TRUE)[2], 
         A = colMeans(param_ests, na.rm = TRUE)[3], B = colMeans(param_ests, na.rm = TRUE)[4], 
         Zc = colMeans(param_ests, na.rm = TRUE)[5]) 

#Compute bias for each parameter
#param_bias_1 <- colMeans(param_ests, na.rm = TRUE) - true_param
param_bias_2 <- colMeans(param_ests, na.rm = TRUE) - true_param


##The results indicate:
# 1. Maximum likelihood estimation of systematic logit GAS is sensitive to the initial parameter selection.
# 1.2. The bias ecountered can probably be significantly reduced by increasing the cros-section and the length of the panel.