#### Purpose: Determine the effect of the panel's length on the calibration of a
#            logit-normal GAS model with a pool of loans, 2 macro variables, one macro factor and one frailty factor


### 2014 - 2016 

#Time-series
tobs <- unique(DM_1C[DateQtr >= 2014]$DateQtr)
#Cross-section
cobs <- matrix(data = 0, nrow = length(tobs), ncol = 1)
for(i in 1:length(tobs)){
  cobs[i] <- nrow(DM_1C[DateQtr == tobs[i]])
}
#Macro 
path_n <- matrix(data = 0, nrow = length(tobs), ncol = 2)
for(i in 1:length(tobs)){
  path_n[i,] <- unique(DM_1C[DateQtr == tobs[i], get(c("HPI", "UR"))])
}
#Default rate
path_l <- matrix(data = 0, nrow = length(tobs), ncol = 1)
for(i in 1:length(tobs)){
  path_l[i] <- sum(DM_1C[DateQtr == tobs[i], Default])/cobs[i]
}

#Likelihood function
loglikelihood <- function(par, path_l, path_n, tobs, cobs){
  f1 <- par[grepl("f1",names(par))]
  Zc <- par[grepl("Zc",names(par))]
  Zm <- cbind(matrix(data = 1, nrow =  ncol(path_n), ncol = 1), matrix(data = 0, nrow =  ncol(path_n), ncol = 1))
  Zm[2,1] <- par[grepl("Zm",names(par))]
  w <- par[grepl("w",names(par))]
  A <- abs(diag(par[grepl("A",names(par))], length(f1), length(f1)))
  B <- diag(par[grepl("B",names(par))], length(f1), length(f1))
  #B <- diag(1/(1+exp(par[grepl("B",names(par))])), length(f1), length(f1))
  ssq <- abs(par[grepl("Sig",names(par))])
  
  Sig <- matrix(data = 0, ncol = ncol(path_n), nrow = ncol(path_n))
  diag(Sig) <- ssq
  Siginv <- diag(1/ssq, ncol = ncol(path_n), nrow = ncol(path_n))
  
  #Common
  score_ <- matrix(data = NA, nrow = length(tobs), ncol = 2)
  f_ <- matrix(data = NA, nrow = length(tobs)+1, ncol = 2)
  loglike <- matrix(data = NA, nrow = length(tobs), ncol = 1)
  #Logit
  p_ <- matrix(data = NA, nrow = length(tobs), ncol = 1)
  score_l_ <- matrix(data = NA, nrow = length(tobs), ncol = 2)
  #Normal
  score_n_ <- matrix(data = NA, nrow = length(tobs), ncol = 2)
  
  #Initialize GAS component
  f_[1,] <- f1
  
  #compute likelihood and other elements at every t
  for(i in 1:length(tobs)){
    #Dynamic probability for logit component
    p_[i] <- 1/(1 + exp(-Zc%*%f_[i,]))
    
    #Score
    score_l_[i,] <- cobs[i]*path_l[i]%*%Zc - cobs[i]*p_[i]%*%Zc
    score_n_[i,] <- Zm%*%Siginv%*%(path_n[i] - Zm%*%f_[i,])
    
    score_[i,] <-  score_l_[i,] + score_n_[i,]
    
    #Log-likelihood
    loglike_l <- cobs[i]*path_l[i]%*%Zc%*%f_[i,] - cobs[i]*log(1 + exp(Zc%*%f_[i,]))
    loglike_n <- -0.5*ncol(path_n)*log(2*pi) - 0.5*log(det(Sig)) - 0.5 * t((path_n[i] - Zm%*%f_[i,])) %*% Siginv %*% (path_n[i] - Zm%*%f_[i,])  
    
    loglike[i] <- loglike_n + loglike_l
    
    f_[i+1,] <- w + A%*%score_[i,] + B%*%f_[i,]
  }
  #compute log-likelihood
  loglike <- sum(loglike)
  return(-loglike)
}

#Initialization
parameters <- c(Zm1 = 0.5, Zc1 = 0.001, Zc2 = 0.001, B1 = 0.5, B2 = 0.5, A1 = 0.5, A2 = 0.5, Sig1 = 1, Sig2 = 1, f1m = 0, 
                f1c = 0) #, w1 = 0, w2 = 0)

#Estimation
fit <- optim(par = parameters, fn = loglikelihood, method = "BFGS" , path_l = path_l, 
             path_n = path_n, cobs = cobs, tobs = tobs, control=list(trace = 1, REPORT=10, maxit = 300), hessian = TRUE)

#Transform parameters to their restricted counter part
par <- c(fit$par[1],  fit$par[2], fit$par[3], fit$par[4], fit$par[5], abs(fit$par[6]), abs(fit$par[7]), abs(fit$par[8]), abs(fit$par[9]),
         fit$par[10], fit$par[11]) #, fit$par[12], fit$par[13])


#Comparison of fitted default rate with observed default rate:
plot(p_, type = "l")
plot(path_l, type = "l")

#Standard Errors
hi <- solve(-fit$hessian)
se <- sqrt(diag(hi))
#P-values
b <- fit$par
zscore <- b / se
p_value <- 2*(1 - pnorm(abs(zscore)))

#Unconditional mean of Factor
Ef <- c(fit$par["w1"], fit$par["w2"])%*%solve((1 - diag(c(fit$par["B1"], fit$par["B2"]) , 2, 2)))

plot(score_[,1], type = "l")
plot(f_[,2], type = "l")




### 2012 - 2016 

#Time-series
tobs <- unique(DM_1C[DateQtr >= 2012]$DateQtr)
#Cross-section
cobs <- matrix(data = 0, nrow = length(tobs), ncol = 1)
for(i in 1:length(tobs)){
  cobs[i] <- nrow(DM_1C[DateQtr == tobs[i]])
}
#Macro 
path_n <- matrix(data = 0, nrow = length(tobs), ncol = 2)
for(i in 1:length(tobs)){
  path_n[i,] <- unique(DM_1C[DateQtr == tobs[i], get(c("HPI", "UR"))])
}
#Default rate
path_l <- matrix(data = 0, nrow = length(tobs), ncol = 1)
for(i in 1:length(tobs)){
  path_l[i] <- sum(DM_1C[DateQtr == tobs[i], Default])/cobs[i]
}


#Estimation
fit <- optim(par = parameters, fn = loglikelihood, method = "BFGS" , path_l = path_l, 
             path_n = path_n, cobs = cobs, tobs = tobs, control=list(trace = 1, REPORT=10, maxit = 300), hessian = TRUE)

#Transform parameters to their restricted counter part
par <- c(fit$par[1],  fit$par[2], fit$par[3], fit$par[4], fit$par[5], abs(fit$par[6]), abs(fit$par[7]), abs(fit$par[8]), abs(fit$par[9]),
         fit$par[10], fit$par[11]) #, fit$par[12], fit$par[13])


#Comparison of fitted default rate with observed default rate:
plot(p_, type = "l")
plot(path_l, type = "l")

#Standard Errors
hi <- solve(-fit$hessian)
se <- sqrt(diag(hi))
#P-values
b <- fit$par
zscore <- b / se
p_value <- 2*(1 - pnorm(abs(zscore)))

#Unconditional mean of Factor
Ef <- c(fit$par["w1"], fit$par["w2"])%*%solve((1 - diag(c(fit$par["B1"], fit$par["B2"]) , 2, 2)))

plot(score_[,2], type = "l")
plot(f_[,2], type = "l")



### 2010 - 2016 

#Time-series
tobs <- unique(DM_1C[DateQtr >= 2010]$DateQtr)
#Cross-section
cobs <- matrix(data = 0, nrow = length(tobs), ncol = 1)
for(i in 1:length(tobs)){
  cobs[i] <- nrow(DM_1C[DateQtr == tobs[i]])
}
#Macro 
path_n <- matrix(data = 0, nrow = length(tobs), ncol = 2)
for(i in 1:length(tobs)){
  path_n[i,] <- unique(DM_1C[DateQtr == tobs[i], get(c("HPI", "UR"))])
}
#Default rate
path_l <- matrix(data = 0, nrow = length(tobs), ncol = 1)
for(i in 1:length(tobs)){
  path_l[i] <- sum(DM_1C[DateQtr == tobs[i], Default])/cobs[i]
}


#Estimation
fit <- optim(par = parameters, fn = loglikelihood, method = "BFGS" , path_l = path_l, 
             path_n = path_n, cobs = cobs, tobs = tobs, control=list(trace = 1, REPORT=10, maxit = 300), hessian = TRUE)

#Transform parameters to their restricted counter part
par <- c(fit$par[1],  fit$par[2], fit$par[3], fit$par[4], fit$par[5], abs(fit$par[6]), abs(fit$par[7]), abs(fit$par[8]), abs(fit$par[9]),
         fit$par[10], fit$par[11]) #, fit$par[12], fit$par[13])


#Comparison of fitted default rate with observed default rate:
plot(p_, type = "l")
plot(path_l, type = "l")

#Standard Errors
hi <- solve(-fit$hessian)
se <- sqrt(diag(hi))
#P-values
b <- fit$par
zscore <- b / se
p_value <- 2*(1 - pnorm(abs(zscore)))

#Unconditional mean of Factor
Ef <- c(fit$par["w1"], fit$par["w2"])%*%solve((1 - diag(c(fit$par["B1"], fit$par["B2"]) , 2, 2)))

plot(score_[,2], type = "l")
plot(f_[,1], type = "l")



### 2008 - 2016 

#Time-series
tobs <- unique(DM_1C[DateQtr >= 2008]$DateQtr)
#Cross-section
cobs <- matrix(data = 0, nrow = length(tobs), ncol = 1)
for(i in 1:length(tobs)){
  cobs[i] <- nrow(DM_1C[DateQtr == tobs[i]])
}
#Macro 
path_n <- matrix(data = 0, nrow = length(tobs), ncol = 2)
for(i in 1:length(tobs)){
  path_n[i,] <- unique(DM_1C[DateQtr == tobs[i], get(c("HPI", "UR"))])
}
#Default rate
path_l <- matrix(data = 0, nrow = length(tobs), ncol = 1)
for(i in 1:length(tobs)){
  path_l[i] <- sum(DM_1C[DateQtr == tobs[i], Default])/cobs[i]
}


parameters <- c(Zm1 = 0.5, Zc1 = 0.001, Zc2 = 0.001, B1 = 0.5, B2 = 0.5, A1 = 0.5, A2 = 0.5, Sig1 = 1, Sig2 = 1, f1m = 0, 
                f1c = 0, w1 = 0, w2 = 0)

#Estimation
fit <- optim(par = parameters, fn = loglikelihood, method = "BFGS" , path_l = path_l, 
             path_n = path_n, cobs = cobs, tobs = tobs, control=list(trace = 1, REPORT=10, maxit = 300), hessian = TRUE)

#Transform parameters to their restricted counter part
par <- c(fit$par[1],  fit$par[2], fit$par[3], fit$par[4], fit$par[4], abs(fit$par[6]), abs(fit$par[7]), abs(fit$par[8]), abs(fit$par[9]),
         fit$par[10], fit$par[11], fit$par[12], fit$par[13])


#Comparison of fitted default rate with observed default rate:
plot(p_, type = "l")
plot(path_l, type = "l")

#Standard Errors
hi <- solve(-fit$hessian)
se <- sqrt(diag(hi))
#P-values
b <- fit$par
zscore <- b / se
p_value <- 2*(1 - pnorm(abs(zscore)))

#Unconditional mean of Factor
Ef <- solve((1 - diag(c(fit$par["B1"], fit$par["B2"]) , 2, 2)))%*%c(fit$par["w1"], fit$par["w2"])

plot(score_[,2], type = "l")
plot(f_[,2], type = "l")


### 2006 - 2016 

#Time-series
tobs <- unique(DM_1C[DateQtr >= 2006]$DateQtr)
#Cross-section
cobs <- matrix(data = 0, nrow = length(tobs), ncol = 1)
for(i in 1:length(tobs)){
  cobs[i] <- nrow(DM_1C[DateQtr == tobs[i]])
}
#Macro 
path_n <- matrix(data = 0, nrow = length(tobs), ncol = 2)
for(i in 1:length(tobs)){
  path_n[i,] <- unique(DM_1C[DateQtr == tobs[i], get(c("HPI", "UR"))])
}
#Default rate
path_l <- matrix(data = 0, nrow = length(tobs), ncol = 1)
for(i in 1:length(tobs)){
  path_l[i] <- sum(DM_1C[DateQtr == tobs[i], Default])/cobs[i]
}


parameters <- c(Zm1 = 0.5, Zc1 = 0.1, Zc2 = 0.001, B1 = 0.3, B2 = 0.2, A1 = 0.5, A2 = 0.5, Sig1 = 0.50, Sig2 = 0.9, f1m = 0, 
                f1c = 0, w1 = 1, w2 = 1)

#Estimation
fit <- optim(par = parameters, fn = loglikelihood, method = "BFGS" , path_l = path_l, 
             path_n = path_n, cobs = cobs, tobs = tobs, control=list(trace = 1, REPORT=10, maxit = 300), hessian = TRUE)

#Transform parameters to their restricted counter part
par <- c(fit$par[1],  fit$par[2], fit$par[3], fit$par[4], fit$par[5], abs(fit$par[6]), abs(fit$par[7]), abs(fit$par[8]), abs(fit$par[9]),
         fit$par[10], fit$par[11], fit$par[12], fit$par[13])


#Comparison of fitted default rate with observed default rate:
plot(p_, type = "l")
plot(path_l, type = "l")

#Standard Errors
hi <- solve(-fit$hessian)
se <- sqrt(diag(hi))
#P-values
b <- fit$par
zscore <- b / se
p_value <- 2*(1 - pnorm(abs(zscore)))

#Unconditional mean of Factor
Ef <- solve((1 - diag(c(fit$par["B1"], fit$par["B2"]) , 2, 2)))%*%c(fit$par["w1"], fit$par["w2"])

plot(score_[,1], type = "l")
plot(f_[,2], type = "l")


### Conclusion:
#1. The logit-normal two factor GAS model completely fails when the length of the panel is short. 
#2. After incorporating more years, the model starts performing well. 
#3. For short periods of time, the factors act in opposite directions. This is an indications that 
#   for short time periods, two factors overfit the model. 
#3. The interpretation of the factors is not straight forward.
#4. The contraint on B has a negative impact on the results. 
#5. After adding the crisis to the panel, the algorithm provided terrible results. 
#5.1 By including the intercept of the factor as a parameter, the results where sensitive again. 
#6. The constraint on B has not been useful in any cases.
#7. The selection of initial parameter values has shown to be of outmost importance to reach a sensitive 
#   maximum of the likelihood funciton. 


