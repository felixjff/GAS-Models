#Purpose: Simulation and estimation of a GAS model with economic cycle. 
#Author: Felix Jose Farias Fueyo


### Simulation ###

simulations <- 1000
tobs <- 120
cobs <- 1000

#GAS process
score <- matrix(data = NA, nrow = tobs, ncol = simulations)
f <- matrix(data = NA, nrow = tobs+1, ncol = simulations)
#Logit
p <- matrix(data = NA, nrow = tobs, ncol = simulations)
gen_success_rate <- matrix(data = NA, nrow = tobs, ncol = simulations)

#Normal - Economic Cycle
y_m <- matrix(data = NA, nrow = tobs ,  ncol = simulations)
sig_nc <- 0.5
mu_nc <- log(12)
sig_c <- 1
mu_c <- log(0.00001)
sig_r <- 2
mu_r <- log(0.1)

# Data Generating Process True Parameter Values
true_param <- c(f1 = 0 , w = 0, alpha = 0.5, beta = 0.5, lambda = 0.3, delta = 0.3)
f1 = 0 ; w = 0; alpha = 0.5; beta = 0.5; lambda = 1; delta = 1

#Level 
tau = 0.008

#simulate
for(i in 1:simulations){
  #initialize factor
  f[1,i] <- true_param["f1"]
  
  y_m[1:36,i]    <- rnorm(n = 36, mean = mu_nc, sd = sig_nc)
  y_m[37:48, i]  <- rnorm(n = 6, mean = mu_c, sd = sig_c)
  y_m[49:56, i]  <- rnorm(n = 8, mean = mu_r, sd = sig_r)
  y_m[57:92, i]  <- rnorm(n = 36, mean = mu_nc, sd = sig_nc)
  y_m[93:104, i] <- rnorm(n = 12, mean = mu_c, sd = sig_c)
  y_m[105:120, i]<- rnorm(n = 16, mean = mu_r, sd = sig_r)
  
  #generate path
  for(t in 1:tobs){
    #Dynamic probability for logit component
    p[t,i] <- tau/(1+ exp(- true_param["w"] - true_param["alpha"]*f[t,i] - true_param["beta"]*y_m[t,i]))
    
    #Instead of storing all observations, only store the generated sucess rate.
    y_l <- rbinom(n = cobs, size = 1, prob = p[t,i])
    gen_success_rate[t,i] <- sum(y_l)/cobs
    
    #scaling comp.
    c <- 1/(1 - tau + exp(- true_param["w"] - true_param["alpha"]*f[t,i] - true_param["beta"]*y_m[t,i]))
    
    #Compute score
    score[t,i] <- (cobs*gen_success_rate[t,i] - cobs*p[t,i]*c)*true_param["alpha"]
    
    f[t+1,i] <- true_param["lambda"]*score[t,i] + true_param["delta"]*f[t,i]
  }
  print(i)
}



### Log-Likelihood ###
loglikelihood <- function(par, path_l, path_n, tobs, cobs){
  f1 <- par[grepl("f1",names(par))]
  w <- par[grepl("w",names(par))]
  lambda <- abs(par[grepl("lambda",names(par))])
  delta <- par[grepl("delta",names(par))]
  alpha <- par[grepl("alpha",names(par))]
  beta <- par[grepl("beta",names(par))]
  
  #Common
  score_ <- matrix(data = NA, nrow = tobs, ncol = 1)
  f_ <- matrix(data = NA, nrow = tobs+1, ncol = 1)
  loglike <- matrix(data = NA, nrow = tobs, ncol = 1)
  #Logit
  p_ <- matrix(data = NA, nrow = tobs, ncol = 1)

  f_[1] <- f1

  for(i in 1:tobs){
    #Estimated PD
    p_[i] <- tau/(1+ exp(- w - alpha*f[i] - beta*path_n[i]))
    
    #scaling comp.
    c <- 1/(1 - tau + exp(- w - alpha*f[i] - beta*y_m[i]))
    
    
    #Score
    score_[i] <- (cobs*path_l[i] - c*cobs*p_[i])*alpha
    
    #Log-likelihood
    loglike[i] <- cobs*path_l[i]*(w + alpha*f[i] + beta*y_m[i]) + log(1 - p_[i])
    
    f_[i+1] <- w + lambda*score_[i] + delta*f_[i]
  }
  loglike <- sum(loglike)
  return(-loglike)
}

#Starting Values
param_0 <- c(f1 = 3 , w = 4, alpha = 0.8, beta = 0.03, lambda = 0.002, delta = 0.001)

#Obtain estimation results for every path
p_e_com_1k <- matrix(data = NA, nrow = simulations, ncol = length(param_0))

for(i in 1:simulations){
  fit <- optim(par = param_0, fn = loglikelihood, method = "BFGS", path_l =  gen_success_rate[,i],
               path_n = y_m[,i], tobs = tobs, cobs = cobs)
  p_e_com_1k[1,] <- fit$par 
  print(i)
}

par <- c(f1 = p_e_com_1k[1,1], w = p_e_com_1k[1,2], alpha = p_e_com_1k[1,3], beta = p_e_com_1k[1,4], lambda = abs(p_e_com_1k[1,5]), delta = p_e_com_1k[1,6])

#Compute bias for each parameter
param_bias <- colMeans(p_e_com_1k, na.rm = TRUE) - true_param


##Conclusion:
#1. Given the choice of the variance and mean of the normal component, the probabilities of default are 
#   mostly (almost only) driven by the normal realizations. Hence, incorporating a GAS component results in a misspecified model 
#   and highly biased maximum likelihood parameters. 