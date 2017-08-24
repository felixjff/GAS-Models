#### Purpose: Determine the effect of the length of the loan panel on the results for the logit and 2 normals GAS model

##### 2014 - 2016 ###

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
  Zm <- par[grepl("Zm",names(par))]
  w <- par[grepl("w",names(par))]
  A <- par[grepl("A",names(par))]
  B <- par[grepl("B",names(par))]
  ssq <- abs(par[grepl("Sig",names(par))])
  
  Sig <- matrix(data = 0, ncol = ncol(path_n), nrow = ncol(path_n))
  diag(Sig) <- ssq
  Siginv <- diag(1/ssq, ncol = ncol(path_n), nrow = ncol(path_n))
  
  #Common
  score_ <- matrix(data = NA, nrow = length(tobs), ncol = 1)
  f_ <- matrix(data = NA, nrow = length(tobs)+1, ncol = 1)
  loglike <- matrix(data = NA, nrow = length(tobs), ncol = 1)
  #Logit
  p_ <- matrix(data = NA, nrow = length(tobs), ncol = 1)
  score_l_ <- matrix(data = NA, nrow = length(tobs), ncol = 1)
  #Normal
  score_n_ <- matrix(data = NA, nrow = length(tobs), ncol = 1)
  
  #Initialize GAS component
  f_[1] <- f1
  
  #compute likelihood and other elements at every t
  for(i in 1:length(tobs)){
    #Dynamic probability for logit component
    p_[i] <- 1/(1 + exp(-Zc*f_[i]))
    
    #Score
    score_l_[i] <- cobs[i]*path_l[i]*Zc - cobs[i]*p_[i]*Zc
    score_n_[i] <- Zm%*%Siginv%*%(path_n[i] - Zm*f_[i])
    
    score_[i] <-  score_l_[i] + score_n_[i]
    
    #Log-likelihood
    loglike_l <- cobs[i]*path_l[i]*Zc*f_[i] - cobs[i]*log(1 + exp(Zc*f_[i]))
    loglike_n <- -0.5*ncol(path_n)*log(2*pi) - 0.5*log(det(Sig)) - 0.5 * t((path_n[i] - Zm*f_[i])) %*% Siginv %*% (path_n[i] - Zm*f_[i])  
    
    loglike[i] <- loglike_n + loglike_l
    
    f_[i+1] <- w + A*score_[i] + B*f_[i]
  }
  #compute log-likelihood
  loglike <- sum(loglike)
  return(-loglike)
}

#Initialization
parameters <- c(Zm1 = 0.5,  Zm2 = 0.5, Zc = 0.001, B = 0.5, A = 0.5, Sig1 = 1, Sig2 = 1, f1 = 0, w = 0)

#Estimation
fit <- optim(par = parameters, fn = loglikelihood, method = "BFGS" , path_l = path_l, 
             path_n = path_n, cobs = cobs, tobs = tobs, control=list(trace = 1, REPORT=1), hessian = TRUE)

#Transform parameters to their restricted counter part
par <- c(fit$par[1],  fit$par[2], fit$par[3], fit$par[4], fit$par[5], abs(fit$par[6]), abs(fit$par[7]), fit$par[8], fit$par[9])

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
Ef <- fit$par["w"]/(1 - fit$par["B"])





##### 2012 - 2016 ###

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
             path_n = path_n, cobs = cobs, tobs = tobs, control=list(trace = 1, REPORT=1), hessian = TRUE)

#Transform parameters to their restricted counter part
par <- c(fit$par[1],  fit$par[2], fit$par[3], fit$par[4], fit$par[5], abs(fit$par[6]), abs(fit$par[7]), fit$par[8], fit$par[9])

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
Ef <- fit$par["w"]/(1 - fit$par["B"])
plot(f_, type = "l")
plot(score_, type = "l")


#Conclusion:
#1. The factor is not centered at its unconditional mean, but the score is centered around 0. 
#2. The generated default rate seems to follow the behaviour observed in the original default rate and be around the same level. 
#3. With small length, the properites of the factor model cannot be directly recognized. 




##### 2010 - 2016 ###

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
             path_n = path_n, cobs = cobs, tobs = tobs, control=list(trace = 1, REPORT=1), hessian = TRUE)

#Transform parameters to their restricted counter part
par <- c(fit$par[1],  fit$par[2], fit$par[3], fit$par[4], fit$par[5], abs(fit$par[6]), abs(fit$par[7]), fit$par[8], fit$par[9])

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
Ef <- fit$par["w"]/(1 - fit$par["B"])
plot(f_, type = "l")
plot(score_, type = "l")


#Conclusion:
#1. The factor is not centered at its unconditional mean, but the score is centered around 0. 
#2. The generated default rate seems to follow the behaviour observed in the original default rate and be around the same level. 
#3. With small length, the properites of the factor model cannot be directly recognized. 




##### 2008 - 2016 ###

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

#Estimation
fit <- optim(par = parameters, fn = loglikelihood, method = "BFGS" , path_l = path_l, 
             path_n = path_n, cobs = cobs, tobs = tobs, control=list(trace = 1, REPORT=1), hessian = TRUE)

#Transform parameters to their restricted counter part
par <- c(fit$par[1],  fit$par[2], fit$par[3], fit$par[4], fit$par[5], abs(fit$par[6]), abs(fit$par[7]), fit$par[8], fit$par[9])

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
Ef <- fit$par["w"]/(1 - fit$par["B"])
plot(f_, type = "l")
plot(score_, type = "l")

#Conclusion:
#1. Factor and score are centered correctly.
#2. The generated default rate seems to follow the behaviour observed in the original default rate and be around the same level. 
#3. Only after a crisis the factor obtains the theoretical properties.  
