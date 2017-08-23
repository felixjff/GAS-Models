#####Purpose: Fititing a systematic-logit GAS model to determine the bias of MLE estimates
###Author: Felix Farias Fueyo


###              ###
###  Simulation  ###
###              ###

#Parameters of true data generating process (DGP)
true_param <- c(f1 = 0.5, w = 0.3, A = 0.5, B = 0.9, Zc = 0.8)
f1 = 0.5; w = 0.3; A = 0.5; B = 0.9; Zc = 0.8

#simulations, cross-section and time-series length
simulations <- 1000
cobs <- 100
tobs <- 100

#Define space for score, factors and realizations for the Binomial indicator variable
score <- matrix(data = NA, nrow = tobs, ncol = simulations)
f <- matrix(data = NA, nrow = tobs, ncol = simulations)
p <- matrix(data = NA, nrow = tobs, ncol = simulations)

#simulate
for(i in 1:simulations){
  #initialize factor
  f[1,i] <- f1
  
  #generate path
  for(t in 1:tobs){
    p[t,i] <- Zc*f[t,i]
    score[t,i] <- p[t,i]*cobs*Zc - cobs*(1/(1+exp(-p[t,i])))*Zc
    f[t+1,i] <- w + A*score[t,i] + B*f[t,i]
  }
  print(i)
}




###              ###
###  Estimation  ###
###              ###

#Define the loglikelihood function for the systematic-logit GAS model 
