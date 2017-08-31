#### Purpose: Calibrate the logit-normal GAS model with a pool of loans increasing length, 2 macro variables, one common latent factor and 
#             5 logit explanatory variables.


#### 2003 - 2016 ###

#Time-series
tobs <- unique(DM_1C[DateQtr > 2003]$DateQtr)
#Cross-section
cobs <- matrix(data = 0, nrow = length(tobs), ncol = 1)
for(i in 1:length(tobs)){
  cobs[i] <- nrow(DM_1C[DateQtr == tobs[i]])
}
#Macro 
path_n <- matrix(data = 0, nrow = length(tobs), ncol = 2)
for(i in 1:length(tobs)){
  path_n[i,1] <- unique(DM_1C[DateQtr == tobs[i], get("HPI")])
  path_n[i,2] <- unique(DM_1C[DateQtr == tobs[i], get("UR")])
}
#Default rate
path_l <- matrix(data = 0, nrow = length(tobs), ncol = 1)
for(i in 1:length(tobs)){
  path_l[i] <- sum(DM_1C[DateQtr == tobs[i], Default])/cobs[i]
}
#Logit Variables
path_l_cscore <- subset(DM_1C, select = c("CSCORE_B","DateQtr"))
path_l_oltv <- subset(DM_1C, select = c("OLTV","DateQtr"))
path_l_dti <- subset(DM_1C, select = c("DTI","DateQtr"))
path_l_origamt <- subset(DM_1C, select = c("ORIG_AMT","DateQtr"))
path_l_lastrt <- subset(DM_1C, select = c("LAST_RT","DateQtr"))


#Likelihood function
loglikelihood <- function(par, path_l, path_n, path_l_cscore, path_l_oltv, path_l_dti, path_l_origamt, path_l_lastrt, tobs, cobs){
  f1 <- par[grepl("f1",names(par))]
  Zc <- par[grepl("Zc",names(par))]
  #Zm <- cbind(matrix(data = 1, nrow =  ncol(path_n), ncol = 1), matrix(data = 0, nrow =  ncol(path_n), ncol = 1))
  #Zm[2,1] <- par[grepl("Zm",names(par))]
  #Zm <- cbind(par[grepl("Zm",names(par))], matrix(data = 0, nrow =  ncol(path_n), ncol = 1))
  Zm <- par[grepl("Zm",names(par))]
  w <- par[grepl("w",names(par))]
  A <- abs(diag(par[grepl("A",names(par))], length(f1), length(f1)))
  B <- diag(par[grepl("B",names(par))], length(f1), length(f1))
  ssq <- abs(par[grepl("Sig",names(par))])
  Beta <- par[grepl("Beta",names(par))]
  
  Sig <- matrix(data = 0, ncol = ncol(path_n), nrow = ncol(path_n))
  diag(Sig) <- ssq
  Siginv <- diag(1/ssq, ncol = ncol(path_n), nrow = ncol(path_n))
  
  #Common
  score_ <- matrix(data = NA, nrow = length(tobs), ncol = 1)
  f_ <- matrix(data = NA, nrow = length(tobs)+1, ncol = 1)
  loglike <- matrix(data = NA, nrow = length(tobs), ncol = 1)
  #Logit
  score_l_ <- matrix(data = NA, nrow = length(tobs), ncol = 1)
  #Normal
  score_n_ <- matrix(data = NA, nrow = length(tobs), ncol = 1)
  
  #Initialize GAS component
  f_[1] <- f1
  
  #compute likelihood and other elements at every t
  for(i in 1:length(tobs)){
    #Dynamic probability for logit component
    sum_cov <- Beta[1]*path_l_cscore[DateQtr == tobs[i]]$CSCORE_B + Beta[2]*path_l_oltv[DateQtr == tobs[i]]$OLTV + 
      Beta[3]*path_l_dti[DateQtr == tobs[i]]$DTI + Beta[4]*path_l_origamt[DateQtr == tobs[i]]$ORIG_AMT +
      Beta[5]*path_l_lastrt[DateQtr == tobs[i]]$LAST_RT
    
    
    p_ <- 1/(1 + exp(Zc*f_[i] - sum_cov))
    
    #Score
    score_l_[i] <- sum(p_)*Zc - cobs[i]*path_l[i]*Zc 
    score_n_[i] <- Zm%*%Siginv%*%(path_n[i,] - Zm*f_[i])
    
    score_[i] <-  score_l_[i] + score_n_[i]
    
    #Log-likelihood
    loglike_l <- -cobs[i]*path_l[i]*Zc*f_[i] - sum(log(1 + exp(sum_cov - Zc*f_[i])))
    loglike_n <- -0.5*ncol(path_n)*log(2*pi) - 0.5*log(det(Sig)) - 0.5 * t((path_n[i,] - Zm*f_[i])) %*% Siginv %*% (path_n[i,] - Zm*f_[i])  
    
    loglike[i] <- loglike_n + loglike_l
    
    f_[i+1] <- w + A*score_[i] + B*f_[i]
  }
  #compute log-likelihood
  loglike <- sum(loglike)
  return(-loglike)
}

#Initialization
parameters <- c(Zm1 = 0.5, Zm2 = 0.5, Zc1 = 0.01, B1 = 0.5, A1 = 0.5, Sig1 = 0.6, Sig2 = 0.6, f1 = 0, Beta = c(1, 1, 1, 1, 1), w1 = 1)

#Estimation
fit <- optim(par = parameters, fn = loglikelihood, method = "BFGS" , path_l = path_l, 
             path_n = path_n, path_l_cscore = path_l_cscore, path_l_oltv = path_l_oltv, path_l_dti = path_l_dti, 
             path_l_origamt = path_l_origamt, path_l_lastrt = path_l_lastrt, 
             cobs = cobs, tobs = tobs, control=list(trace = 1, REPORT=10, maxit = 300))

#Transform parameters to their restricted counter part
par <- c(fit$par[1],  fit$par[2], fit$par[3], fit$par[4], fit$par[5], abs(fit$par[6]), abs(fit$par[7]), abs(fit$par[8]), abs(fit$par[9]),
         fit$par[10], fit$par[11], fit$par[12], fit$par[13], fit$par[14])

par_03_16 <- par

#Comparison of fitted default rate with observed default rate:
Beta_ <- c(fit$par["Beta1"], fit$par["Beta2"], fit$par["Beta3"], fit$par["Beta4"], fit$par["Beta5"])
p_average <- matrix(data = NA, nrow = length(tobs), ncol = 1)
for(i in 1:length(tobs)){
  #Dynamic probability for logit component
  sum_cov <- Beta_[1]*path_l_cscore[DateQtr == tobs[i]]$CSCORE_B + Beta_[2]*path_l_oltv[DateQtr == tobs[i]]$OLTV + 
    Beta_[3]*path_l_dti[DateQtr == tobs[i]]$DTI + Beta_[4]*path_l_origamt[DateQtr == tobs[i]]$ORIG_AMT +
    Beta_[5]*path_l_lastrt[DateQtr == tobs[i]]$LAST_RT
  
  
  p_average[i] <- sum(1/(1 + exp(fit$par["Zc1"]*f_[i] - sum_cov)))/cobs[i]
}

plot(p_average, type = "l")
plot(path_l, type = "l")

#Unconditional mean of Factor
Ef <- fit$par["w1"]/(1 - fit$par["B1"])

plot(score_, type = "l")
plot(f_, type = "l")
f_03_16 <- f_
score_03_16 <- score_

#Standard Errors
hi <- numDeriv::hessian(x = par, func = loglikelihood, method= "Richardson", path_l = path_l, 
                        path_n = path_n, path_l_cscore = path_l_cscore, path_l_oltv = path_l_oltv, path_l_dti = path_l_dti, 
                        path_l_origamt = path_l_origamt, path_l_lastrt = path_l_lastrt, 
                        cobs = cobs, tobs = tobs)

hi <- solve(-hi)
se <- sqrt(diag(hi))
se_03_16 <- se
#P-values
b <- fit$par
zscore <- b / se
p_value <- 2*(1 - pnorm(abs(zscore)))
p_value_03_16 <- p_value



#### 2006 - 2016 ###

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
  path_n[i,1] <- unique(DM_1C[DateQtr == tobs[i], get("HPI")])
  path_n[i,2] <- unique(DM_1C[DateQtr == tobs[i], get("UR")])
}
#Default rate
path_l <- matrix(data = 0, nrow = length(tobs), ncol = 1)
for(i in 1:length(tobs)){
  path_l[i] <- sum(DM_1C[DateQtr == tobs[i], Default])/cobs[i]
}

#Initialization
parameters <- c(Zm1 = 0.5, Zm2 = 0.5, Zc1 = 0.01, B1 = 0.5, A1 = 0.5, Sig1 = 0.6, Sig2 = 0.6, f1 = 0, Beta = c(1, 1, 1, 1, 1), w1 = 1)

#Estimation
fit <- optim(par = parameters, fn = loglikelihood, method = "BFGS" , path_l = path_l, 
             path_n = path_n, path_l_cscore = path_l_cscore, path_l_oltv = path_l_oltv, path_l_dti = path_l_dti, 
             path_l_origamt = path_l_origamt, path_l_lastrt = path_l_lastrt, 
             cobs = cobs, tobs = tobs, control=list(trace = 1, REPORT=10, maxit = 300))

#Transform parameters to their restricted counter part
par <- c(fit$par[1],  fit$par[2], fit$par[3], fit$par[4], fit$par[5], abs(fit$par[6]), abs(fit$par[7]), abs(fit$par[8]), abs(fit$par[9]),
         fit$par[10], fit$par[11], fit$par[12], fit$par[13], fit$par[14])

par_06_16 <- par

#Comparison of fitted default rate with observed default rate:
Beta_ <- c(fit$par["Beta1"], fit$par["Beta2"], fit$par["Beta3"], fit$par["Beta4"], fit$par["Beta5"])
p_average <- matrix(data = NA, nrow = length(tobs), ncol = 1)
for(i in 1:length(tobs)){
  #Dynamic probability for logit component
  sum_cov <- Beta_[1]*path_l_cscore[DateQtr == tobs[i]]$CSCORE_B + Beta_[2]*path_l_oltv[DateQtr == tobs[i]]$OLTV + 
    Beta_[3]*path_l_dti[DateQtr == tobs[i]]$DTI + Beta_[4]*path_l_origamt[DateQtr == tobs[i]]$ORIG_AMT +
    Beta_[5]*path_l_lastrt[DateQtr == tobs[i]]$LAST_RT
  
  
  p_average[i] <- sum(1/(1 + exp(fit$par["Zc1"]*f_[i] - sum_cov)))/cobs[i]
}


plot(p_average, type = "l")
plot(path_l, type = "l")

#Unconditional mean of Factor
Ef <- fit$par["w1"]/(1 - fit$par["B1"])

plot(score_, type = "l")
plot(f_, type = "l")
f_06_16 <- f_
score_06_16 <- score_

#Standard Errors
hi <- solve(-fit$hessian)
se <- sqrt(diag(hi))
se_06_16 <- se
#P-values
b <- fit$par
zscore <- b / se
p_value <- 2*(1 - pnorm(abs(zscore)))
p_value_06_16 <- p_value




#### 2008 - 2016 ###

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
  path_n[i,1] <- unique(DM_1C[DateQtr == tobs[i], get("HPI")])
  path_n[i,2] <- unique(DM_1C[DateQtr == tobs[i], get("UR")])
}
#Default rate
path_l <- matrix(data = 0, nrow = length(tobs), ncol = 1)
for(i in 1:length(tobs)){
  path_l[i] <- sum(DM_1C[DateQtr == tobs[i], Default])/cobs[i]
}

#Initialization
parameters <- c(Zm1 = 0.5, Zm2 = 0.5, Zc1 = 0.01, B1 = 0.5, A1 = 0.5, Sig1 = 0.6, Sig2 = 0.6, f1 = 0, Beta = c(1, 1, 1, 1, 1), w1 = 1)

#Estimation
fit <- optim(par = parameters, fn = loglikelihood, method = "BFGS" , path_l = path_l, 
             path_n = path_n, path_l_cscore = path_l_cscore, path_l_oltv = path_l_oltv, path_l_dti = path_l_dti, 
             path_l_origamt = path_l_origamt, path_l_lastrt = path_l_lastrt, 
             cobs = cobs, tobs = tobs, control=list(trace = 1, REPORT=10, maxit = 300))

#Transform parameters to their restricted counter part
par <- c(fit$par[1],  fit$par[2], fit$par[3], fit$par[4], fit$par[5], abs(fit$par[6]), abs(fit$par[7]), abs(fit$par[8]), abs(fit$par[9]),
         fit$par[10], fit$par[11], fit$par[12], fit$par[13], fit$par[14])

par_08_16 <- par

#Comparison of fitted default rate with observed default rate:
Beta_ <- c(fit$par["Beta1"], fit$par["Beta2"], fit$par["Beta3"], fit$par["Beta4"], fit$par["Beta5"])
p_average <- matrix(data = NA, nrow = length(tobs), ncol = 1)
for(i in 1:length(tobs)){
  #Dynamic probability for logit component
  sum_cov <- Beta_[1]*path_l_cscore[DateQtr == tobs[i]]$CSCORE_B + Beta_[2]*path_l_oltv[DateQtr == tobs[i]]$OLTV + 
    Beta_[3]*path_l_dti[DateQtr == tobs[i]]$DTI + Beta_[4]*path_l_origamt[DateQtr == tobs[i]]$ORIG_AMT +
    Beta_[5]*path_l_lastrt[DateQtr == tobs[i]]$LAST_RT
  
  
  p_average[i] <- sum(1/(1 + exp(fit$par["Zc1"]*f_[i] - sum_cov)))/cobs[i]
}


plot(p_average, type = "l")
plot(path_l, type = "l")

#Unconditional mean of Factor
Ef <- fit$par["w1"]/(1 - fit$par["B1"])

plot(score_, type = "l")
plot(f_, type = "l")
f_08_16 <- f_
score_08_16 <- score_

#Standard Errors
hi <- solve(-fit$hessian)
se <- sqrt(diag(hi))
se_08_16 <- se
#P-values
b <- fit$par
zscore <- b / se
p_value <- 2*(1 - pnorm(abs(zscore)))
p_value_08_16 <- p_value



#### 2010 - 2016 ###

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
  path_n[i,1] <- unique(DM_1C[DateQtr == tobs[i], get("HPI")])
  path_n[i,2] <- unique(DM_1C[DateQtr == tobs[i], get("UR")])
}
#Default rate
path_l <- matrix(data = 0, nrow = length(tobs), ncol = 1)
for(i in 1:length(tobs)){
  path_l[i] <- sum(DM_1C[DateQtr == tobs[i], Default])/cobs[i]
}

#Initialization
parameters <- c(Zm1 = 0.5, Zm2 = 0.5, Zc1 = 0.01, B1 = 0.5, A1 = 0.5, Sig1 = 0.6, Sig2 = 0.6, f1 = 0, Beta = c(1, 1, 1, 1, 1), w1 = 1)

#Estimation
fit <- optim(par = parameters, fn = loglikelihood, method = "BFGS" , path_l = path_l, 
             path_n = path_n, path_l_cscore = path_l_cscore, path_l_oltv = path_l_oltv, path_l_dti = path_l_dti, 
             path_l_origamt = path_l_origamt, path_l_lastrt = path_l_lastrt, 
             cobs = cobs, tobs = tobs, control=list(trace = 1, REPORT=10, maxit = 300))

#Transform parameters to their restricted counter part
par <- c(fit$par[1],  fit$par[2], fit$par[3], fit$par[4], fit$par[5], abs(fit$par[6]), abs(fit$par[7]), abs(fit$par[8]), abs(fit$par[9]),
         fit$par[10], fit$par[11], fit$par[12], fit$par[13], fit$par[14])

par_10_16 <- par

#Comparison of fitted default rate with observed default rate:
Beta_ <- c(fit$par["Beta1"], fit$par["Beta2"], fit$par["Beta3"], fit$par["Beta4"], fit$par["Beta5"])
p_average <- matrix(data = NA, nrow = length(tobs), ncol = 1)
for(i in 1:length(tobs)){
  #Dynamic probability for logit component
  sum_cov <- Beta_[1]*path_l_cscore[DateQtr == tobs[i]]$CSCORE_B + Beta_[2]*path_l_oltv[DateQtr == tobs[i]]$OLTV + 
    Beta_[3]*path_l_dti[DateQtr == tobs[i]]$DTI + Beta_[4]*path_l_origamt[DateQtr == tobs[i]]$ORIG_AMT +
    Beta_[5]*path_l_lastrt[DateQtr == tobs[i]]$LAST_RT
  
  
  p_average[i] <- sum(1/(1 + exp(fit$par["Zc1"]*f_[i] - sum_cov)))/cobs[i]
}

plot(p_average, type = "l")
plot(path_l, type = "l")

#Unconditional mean of Factor
Ef <- fit$par["w1"]/(1 - fit$par["B1"])

plot(score_, type = "l")
plot(f_, type = "l")
f_10_16 <- f_
score_10_16 <- score_

#Standard Errors
hi <- solve(-fit$hessian)
se <- sqrt(diag(hi))
se_10_16 <- se
#P-values
b <- fit$par
zscore <- b / se
p_value <- 2*(1 - pnorm(abs(zscore)))
p_value_10_16 <- p_value




#### 2012 - 2016 ###

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
  path_n[i,1] <- unique(DM_1C[DateQtr == tobs[i], get("HPI")])
  path_n[i,2] <- unique(DM_1C[DateQtr == tobs[i], get("UR")])
}
#Default rate
path_l <- matrix(data = 0, nrow = length(tobs), ncol = 1)
for(i in 1:length(tobs)){
  path_l[i] <- sum(DM_1C[DateQtr == tobs[i], Default])/cobs[i]
}

#Initialization
parameters <- c(Zm1 = 0.5, Zm2 = 0.5, Zc1 = 0.01, B1 = 0.5, A1 = 0.5, Sig1 = 0.6, Sig2 = 0.6, f1 = 0, Beta = c(1, 1, 1, 1, 1), w1 = 1)

#Estimation
fit <- optim(par = parameters, fn = loglikelihood, method = "BFGS" , path_l = path_l, 
             path_n = path_n, path_l_cscore = path_l_cscore, path_l_oltv = path_l_oltv, path_l_dti = path_l_dti, 
             path_l_origamt = path_l_origamt, path_l_lastrt = path_l_lastrt, 
             cobs = cobs, tobs = tobs, control=list(trace = 1, REPORT=10, maxit = 300))

#Transform parameters to their restricted counter part
par <- c(fit$par[1],  fit$par[2], fit$par[3], fit$par[4], fit$par[5], abs(fit$par[6]), abs(fit$par[7]), abs(fit$par[8]), abs(fit$par[9]),
         fit$par[10], fit$par[11], fit$par[12], fit$par[13], fit$par[14])

par_12_16 <- par

#Comparison of fitted default rate with observed default rate:
Beta_ <- c(fit$par["Beta1"], fit$par["Beta2"], fit$par["Beta3"], fit$par["Beta4"], fit$par["Beta5"])
p_average <- matrix(data = NA, nrow = length(tobs), ncol = 1)
for(i in 1:length(tobs)){
  #Dynamic probability for logit component
  sum_cov <- Beta_[1]*path_l_cscore[DateQtr == tobs[i]]$CSCORE_B + Beta_[2]*path_l_oltv[DateQtr == tobs[i]]$OLTV + 
    Beta_[3]*path_l_dti[DateQtr == tobs[i]]$DTI + Beta_[4]*path_l_origamt[DateQtr == tobs[i]]$ORIG_AMT +
    Beta_[5]*path_l_lastrt[DateQtr == tobs[i]]$LAST_RT
  
  
  p_average[i] <- sum(1/(1 + exp(fit$par["Zc1"]*f_[i] - sum_cov)))/cobs[i]
}

plot(p_average, type = "l")
plot(path_l, type = "l")

#Unconditional mean of Factor
Ef <- fit$par["w1"]/(1 - fit$par["B1"])

plot(score_, type = "l")
plot(f_, type = "l")
f_12_16 <- f_
score_12_16 <- score_

#Standard Errors
hi <- solve(-fit$hessian)
se <- sqrt(diag(hi))
se_12_16 <- se
#P-values
b <- fit$par
zscore <- b / se
p_value <- 2*(1 - pnorm(abs(zscore)))
p_value_12_16 <- p_value




#### 2014 - 2016 ###

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
  path_n[i,1] <- unique(DM_1C[DateQtr == tobs[i], get("HPI")])
  path_n[i,2] <- unique(DM_1C[DateQtr == tobs[i], get("UR")])
}
#Default rate
path_l <- matrix(data = 0, nrow = length(tobs), ncol = 1)
for(i in 1:length(tobs)){
  path_l[i] <- sum(DM_1C[DateQtr == tobs[i], Default])/cobs[i]
}

#Initialization
parameters <- c(Zm1 = 0.5, Zm2 = 0.5, Zc1 = 0.01, B1 = 0.5, A1 = 0.5, Sig1 = 0.6, Sig2 = 0.6, f1 = 0, Beta = c(1, 1, 1, 1, 1), w1 = 1)

#Estimation
fit <- optim(par = parameters, fn = loglikelihood, method = "BFGS" , path_l = path_l, 
             path_n = path_n, path_l_cscore = path_l_cscore, path_l_oltv = path_l_oltv, path_l_dti = path_l_dti, 
             path_l_origamt = path_l_origamt, path_l_lastrt = path_l_lastrt, 
             cobs = cobs, tobs = tobs, control=list(trace = 1, REPORT=10, maxit = 300))

#Transform parameters to their restricted counter part
par <- c(fit$par[1],  fit$par[2], fit$par[3], fit$par[4], fit$par[5], abs(fit$par[6]), abs(fit$par[7]), abs(fit$par[8]), abs(fit$par[9]),
         fit$par[10], fit$par[11], fit$par[12], fit$par[13], fit$par[14])

par_14_16 <- par

#Comparison of fitted default rate with observed default rate:
Beta_ <- c(fit$par["Beta1"], fit$par["Beta2"], fit$par["Beta3"], fit$par["Beta4"], fit$par["Beta5"])
p_average <- matrix(data = NA, nrow = length(tobs), ncol = 1)
for(i in 1:length(tobs)){
  #Dynamic probability for logit component
  sum_cov <- Beta_[1]*path_l_cscore[DateQtr == tobs[i]]$CSCORE_B + Beta_[2]*path_l_oltv[DateQtr == tobs[i]]$OLTV + 
    Beta_[3]*path_l_dti[DateQtr == tobs[i]]$DTI + Beta_[4]*path_l_origamt[DateQtr == tobs[i]]$ORIG_AMT +
    Beta_[5]*path_l_lastrt[DateQtr == tobs[i]]$LAST_RT
  
  
  p_average[i] <- sum(1/(1 + exp(fit$par["Zc1"]*f_[i] - sum_cov)))/cobs[i]
}

plot(p_average, type = "l")
plot(path_l, type = "l")

#Unconditional mean of Factor
Ef <- fit$par["w1"]/(1 - fit$par["B1"])

plot(score_, type = "l")
plot(f_, type = "l")
f_14_16 <- f_
score_14_16 <- score_

#Standard Errors
hi <- solve(-fit$hessian)
se <- sqrt(diag(hi))
se_14_16 <- se
#P-values
b <- fit$par
zscore <- b / se
p_value <- 2*(1 - pnorm(abs(zscore)))
p_value_14_16 <- p_value