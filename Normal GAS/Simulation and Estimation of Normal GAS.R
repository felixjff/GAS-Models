#####Purpose: Fititing a normal GAS with a normal GAS model to determine the bias of MLE estimates
###Author: Felix Farias Fueyo


###              ###
###  Simulation  ###
###              ###

###Generate the normal GAS sample
simulations <- 1000
tobs <- 100
 
# Data Generating Process True Parameter Values
true_param <- c(f1 = 0 , Zm = 0.5, w = 0, A = 0.5, B = 0.5, Sigma_squared = 1)
f1 = 0 ; Zm = 0.5; w = 0; A = 0.5; B = 0.5; Sigma_squared = 1

#Storage value for scores, factors and normal variable
y <- matrix(data = NA, nrow = tobs ,  ncol = simulations)
score <- matrix(data = NA, nrow = tobs, ncol = simulations)
f <- matrix(data = NA, nrow = tobs + 1, ncol = simulations)

for(i in 1:simulations){
 #initialize factor
  f[1,i] <- f1
  
  #generate path
  for(t in 1:tobs){
    y[t,i] <- rnorm(1, mean = Zm*f[t,i], sd = Sigma_squared)
    score[t,i] <- (1/Sigma_squared)*Zm*(y[t,i] - Zm*f[t,i])
    f[t+1,i] <- w + A*score[t,i] + B*f[t,i]
  }
  print(i)
}



###              ###
###  Estimation  ###
###              ###

#Define the loglikelihood function for the normal gas model
loglikelihood <- function(param, path, tobs){
    f1 <- param[grepl("f1",names(param))]
    Zm <- param[grepl("Zm",names(param))]
    w <- param[grepl("w",names(param))]
    A <- param[grepl("A",names(param))]
    B <- param[grepl("B",names(param))]
    ssq <- param[grepl("Sigma_squared",names(param))]
    
    # create storage space for the score, GAS component and the log likelihood
    score <- matrix(data = NA, nrow = tobs, ncol = 1)
    loglike <- matrix(data = NA, nrow = tobs, ncol = 1)
    f <- matrix(data = NA, nrow = tobs + 1, ncol = 1)
    
    #Initialize GAS component
    f[1] <- f1
    
    #compute likelihood and other elements at every t
    for(i in 1:tobs){
      score[i] <- (1/ssq)*Zm*(path[i] - Zm*f[i])
      loglike[i] <- -0.5*log(2*pi*ssq) - 0.5*(1/ssq)*(path[i] - Zm*f[i])^2
      f[i+1] <- w + A*score[i] + B*f[i]
    }
    
    #compute log-likelihood
    loglike <- sum(loglike)
    return(-loglike)
}


#Obtain estimation results for every path
param_0 <- c(f1 = 1 , Zm = 2, w = 1, A = 0.1, B = 0.99, Sigma_squared = 3)

param_ests <- matrix(data = NA, nrow = simulations, ncol = length(param_0))

for(i in 1:simulations){
  fit <- optim(par = param_0, fn = loglikelihood, , method = "BFGS", path = y[,i], tobs = tobs)
  param_ests[i,] <- fit$par 
  print(i)
}


#Compute bias for each parameter
param_bias <- colMeans(param_ests) - true_param
