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

#Create storage space for the score, factor, normal variable, success indicator, true success rate and generated success rate
#Common
score_m <- matrix(data = NA, nrow = tobs, ncol = simulations)
score_c <- matrix(data = NA, nrow = tobs, ncol = simulations)
f <- matrix(data = NA, nrow = tobs+1, ncol = simulations)
#Logit
p <- matrix(data = NA, nrow = tobs, ncol = simulations)
gen_success_rate <- matrix(data = NA, nrow = tobs, ncol = simulations)
#Normal
y <- matrix(data = NA, nrow = tobs ,  ncol = simulations)

# Data Generating Process True Parameter Values
true_param <- c(f1 = 0 , Zm = 0.5, Zc = 0.5, w = 0, A = 0.5, B = 0.5, Sigma_squared = 1)
f1 = 0 ; Zm = 0.5; w = 0; A = 0.5; B = 0.5; Sigma_squared = 1; Zc = 0.5

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