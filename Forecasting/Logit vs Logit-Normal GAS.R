#### Purpose: Compare the performance of a simple logit model with observable macro risk factor 
#    with the perfromance of logit-normal one factor GAS model. The models are used to forecast a crisis and a stable grwoth period.



#### LOGIT MODEL ###

## Log-Likelihood
#Step2
Loglikelihood_L1 <- function(beta, data, qrts){
  logl <- 0
  
  for(i in qrts){
    temp <- as.data.table(data[DateQtr == i])
    temp[Default == 1, Lcont := log(exp(beta[1] + beta[2]*CSCORE_B + beta[3]*OLTV + beta[4]*DTI +
                                          beta[5]*ORIG_AMT + beta[6]*LAST_RT)/(1+exp(beta[1] + beta[2]*CSCORE_B + beta[3]*OLTV + beta[4]*DTI +
                                                                                       beta[5]*ORIG_AMT + beta[6]*LAST_RT)))]
    temp[Default == 0, Lcont := log(1 - exp(beta[1] + beta[2]*CSCORE_B + beta[3]*OLTV + beta[4]*DTI +
                                              beta[5]*ORIG_AMT + beta[6]*LAST_RT)/(1+exp(beta[1] + beta[2]*CSCORE_B + beta[3]*OLTV + beta[4]*DTI +
                                                                                           beta[5]*ORIG_AMT + beta[6]*LAST_RT)))]
    logl <- logl + sum(temp$Lcont)
  }
  
  return(-logl)
}

#Step 1
Loglikelihood_L2 <- function(par_pit, data, qrts){
  logl <- 0
  
  for(i in qrts){
    temp <- as.data.table(data[DateQtr == i])
    temp[Default == 1, Lcont := log(exp(par_pit[1]*p_it_star + par_pit[2]*GDP + par_pit[3]*HPI + par_pit[4]*UR + par_pit[5]*IPG)
                                    /(1+exp(par_pit[1]*p_it_star + par_pit[2]*GDP + par_pit[3]*HPI + par_pit[4]*UR + par_pit[5]*IPG)))]
    temp[Default == 0, Lcont := log(1 - exp(par_pit[1]*p_it_star + par_pit[2]*GDP + par_pit[3]*HPI + par_pit[4]*UR + par_pit[5]*IPG)
                                    /(1+exp(par_pit[1]*p_it_star + par_pit[2]*GDP + par_pit[3]*HPI + par_pit[4]*UR 
                                            + par_pit[5]*IPG)))]
    logl <- logl + sum(temp$Lcont)
  }
  
  return(-logl)
}


## 2003 - 2007
#Calibration: Step 1
beta <- matrix(1, 1, 6)
tobs <- unique(DM_1C[DateQtr >= 2003 & DateQtr < 2008, DateQtr])

fit <- optim(par = beta, fn = Loglikelihood_L1, method = "BFGS" , qrts = tobs, data = DM_1C[DateQtr >= 2003 & DateQtr < 2008],  
                   control=list(trace=TRUE), hessian = TRUE)

hi <- solve(-fit$hessian)
se <- sqrt(diag(hi))
se_03_07_S1 <- se

par <- fit$par
par_03_07_S1 <- par

b <- par
zscore <- b / se
p_value <- 2*(1 - pnorm(abs(zscore)))
pvalue_03_07_S1 <- p_value

p_it_star <- b[1] + b[2]*DM_1C[DateQtr >= 2003 & DateQtr < 2008]$CSCORE_B + b[3]*DM_1C[DateQtr >= 2003 & DateQtr < 2008]$OLTV + 
             b[4]*DM_1C[DateQtr >= 2003 & DateQtr < 2008]$DTI + b[5]*DM_1C[DateQtr >= 2003 & DateQtr < 2008]$ORIG_AMT + 
             b[6]*DM_1C[DateQtr >= 2003 & DateQtr < 2008]$LAST_RT 

DM_1C[, p_it_star := NA_real_]
DM_1C[DateQtr >= 2003 & DateQtr < 2008, p_it_star := b[1] + b[2]*CSCORE_B + b[3]*OLTV + b[4]*DTI + b[5]*ORIG_AMT + b[6]*LAST_RT]

#Calibration: Step 2
fit <- optim(par = beta, fn = Loglikelihood_L2, method = "BFGS" , qrts = tobs, data = DM_1C[DateQtr >= 2003 & DateQtr < 2008],  
                    control=list(trace=TRUE), hessian = TRUE)

hi<- solve(-fit$hessian)
se <- sqrt(diag(hi))
se_03_07_S2 <- se

par <- fit$par
par_03_07_S2 <- par

b <- par 
zscore <- b / se
p_value <- 2*(1 - pnorm(abs(zscore)))
pvalue_03_07_S1 <- p_value

#Forecast
#***Use only param. at 1% sig. level
DM_1C[, PD_forecast_L1 := NA_real_]
DM_1C[DateQtr >= 2008 & DateQtr < 2010, p_it_star := par_03_07_S1[1] + par_03_07_S1[2]*CSCORE_B + par_03_07_S1[6]*LAST_RT]
DM_1C[DateQtr >= 2008 & DateQtr < 2010, PD_forecast_L := exp(par_03_07_S2[1]*p_it_star)/(1 + exp(par_03_07_S2[1]*p_it_star))]

# Performance
#AUC
auc_03_07_L_1 <- auc(DM_1C[DateQtr == 2008.00]$Default, DM_1C[DateQtr == 2008.00]$PD_forecast_L)
auc_03_07_L_2 <- auc(DM_1C[DateQtr == 2008.25]$Default, DM_1C[DateQtr == 2008.25]$PD_forecast_L)
auc_03_07_L_3 <- auc(DM_1C[DateQtr == 2008.50]$Default, DM_1C[DateQtr == 2008.50]$PD_forecast_L)
auc_03_07_L_4 <- auc(DM_1C[DateQtr == 2008.75]$Default, DM_1C[DateQtr == 2008.75]$PD_forecast_L)
auc_03_07_L_5 <- auc(DM_1C[DateQtr == 2009.00]$Default, DM_1C[DateQtr == 2009.00]$PD_forecast_L)
auc_03_07_L_6 <- auc(DM_1C[DateQtr == 2009.25]$Default, DM_1C[DateQtr == 2009.25]$PD_forecast_L)
auc_03_07_L_7 <- auc(DM_1C[DateQtr == 2009.50]$Default, DM_1C[DateQtr == 2009.50]$PD_forecast_L)
auc_03_07_L_8 <- auc(DM_1C[DateQtr == 2009.75]$Default, DM_1C[DateQtr == 2009.75]$PD_forecast_L)

#MAD
tobs <- unique(DM_1C[DateQtr >= 2008 & DateQtr < 2010]$DateQtr)
def_rate <- matrix(data = NA, nrow = length(tobs), ncol = 1)
def_rate_estim <- matrix(data = NA, nrow = length(tobs), ncol = 1)

for(t in 1:length(tobs)){
  def_rate_estim[t] <- sum(DM_1C[DateQtr == tobs[t]]$PD_forecast_L)/nrow(DM_1C[DateQtr == tobs[t]])
  def_rate[t] <- sum(DM_1C[DateQtr == tobs[t]]$Default)/nrow(DM_1C[DateQtr == tobs[t]])
}

MAD_03_07_L <- sum(abs(def_rate - def_rate_estim))/length(tobs)

DM_1C[, p_it_star := NA_real_]
DM_1C[, PD_forecast_L := NA_real_]


## 2003 - 2014
#Calibration: Step 1
beta <- matrix(1, 1, 6)
tobs <- unique(DM_1C[DateQtr >= 2003 & DateQtr < 2014, DateQtr])

fit <- optim(par = beta, fn = Loglikelihood_L1, method = "BFGS" , qrts = tobs, data = DM_1C[DateQtr >= 2003 & DateQtr < 2014],  
             control=list(trace=TRUE), hessian = TRUE)

hi <- solve(-fit$hessian)
se <- sqrt(diag(hi))
se_03_14_S2 <- se

par <- fit$par
par_03_14_S1 <- par

b <- par
zscore <- b / se
p_value <- 2*(1 - pnorm(abs(zscore)))
pvalue_03_14_S1 <- p_value

p_it_star <- b[1] + b[2]*DM_1C[DateQtr >= 2003 & DateQtr < 2014]$CSCORE_B + b[3]*DM_1C[DateQtr >= 2003 & DateQtr < 2014]$OLTV + 
  b[4]*DM_1C[DateQtr >= 2003 & DateQtr < 2014]$DTI + b[5]*DM_1C[DateQtr >= 2003 & DateQtr < 2014]$ORIG_AMT + 
  b[6]*DM_1C[DateQtr >= 2003 & DateQtr < 2014]$LAST_RT 

DM_1C[, p_it_star := NA_real_]
DM_1C[DateQtr >= 2003 & DateQtr < 2014, p_it_star := b[1] + b[2]*CSCORE_B + b[3]*OLTV + b[4]*DTI + b[5]*ORIG_AMT + b[6]*LAST_RT]

#Calibration: Step 2
fit <- optim(par = beta, fn = Loglikelihood_L2, method = "BFGS" , qrts = tobs, data = DM_1C[DateQtr >= 2003 & DateQtr < 2014],  
                    control=list(trace=TRUE), hessian = TRUE)

hi<- solve(-fit$hessian)
se <- sqrt(diag(hi))
se_03_14_S2 <- se

par <- fit$par
par_03_14_S2 <- par

b <- par 
zscore <- b / se
p_value <- 2*(1 - pnorm(abs(zscore)))
pvalue_03_14_S1 <- p_value


#Forecast
#***Use only param. at 1% sig. level
DM_1C[, PD_forecast_L := NA_real_]
DM_1C[DateQtr >= 2014 & DateQtr < 2016, p_it_star := par_03_14_S1[1] + par_03_14_S1[2]*CSCORE_B + 
        #par_03_14_S1[3]*OLTV + par_03_14_S1[4]*DTI + par_03_14_S1[2]*ORIG_AMT + 
        par_03_14_S1[6]*LAST_RT]
DM_1C[DateQtr >= 2014 & DateQtr < 2016, PD_forecast_L := exp(par_03_07_S2[1]*p_it_star + par_03_07_S2[4]*UR)/
        (1 + exp(par_03_07_S2[1]*p_it_star + par_03_07_S2[4]*UR))]

# Performance
#AUC
auc_14_16_L_1 <- auc(DM_1C[DateQtr == 2014.00]$Default, DM_1C[DateQtr == 2014.00]$PD_forecast_L)
auc_14_16_L_2 <- auc(DM_1C[DateQtr == 2014.25]$Default, DM_1C[DateQtr == 2014.25]$PD_forecast_L)
auc_14_16_L_3 <- auc(DM_1C[DateQtr == 2014.50]$Default, DM_1C[DateQtr == 2014.50]$PD_forecast_L)
auc_14_16_L_4 <- auc(DM_1C[DateQtr == 2014.75]$Default, DM_1C[DateQtr == 2014.75]$PD_forecast_L)
auc_14_16_L_5 <- auc(DM_1C[DateQtr == 2015.00]$Default, DM_1C[DateQtr == 2015.00]$PD_forecast_L)
auc_14_16_L_6 <- auc(DM_1C[DateQtr == 2015.25]$Default, DM_1C[DateQtr == 2015.25]$PD_forecast_L)
auc_14_16_L_7 <- auc(DM_1C[DateQtr == 2015.50]$Default, DM_1C[DateQtr == 2015.50]$PD_forecast_L)
auc_14_16_L_8 <- auc(DM_1C[DateQtr == 2015.75]$Default, DM_1C[DateQtr == 2015.75]$PD_forecast_L)

#MAD
tobs <- unique(DM_1C[DateQtr >= 2014 & DateQtr < 2016]$DateQtr)
def_rate <- matrix(data = NA, nrow = length(tobs), ncol = 1)
def_rate_estim <- matrix(data = NA, nrow = length(tobs), ncol = 1)

for(t in 1:length(tobs)){
  def_rate_estim[t] <- sum(DM_1C[DateQtr == tobs[t]]$PD_forecast_L)/nrow(DM_1C[DateQtr == tobs[t]])
  def_rate[t] <- sum(DM_1C[DateQtr == tobs[t]]$Default)/nrow(DM_1C[DateQtr == tobs[t]])
}

MAD_03_14_L <- sum(abs(def_rate - def_rate_estim))/length(tobs)

DM_1C[, p_it_star := NA_real_]
DM_1C[, PD_forecast_L := NA_real_]



#### GAS MODEL ###

## Log-Likelihood
#Likelihood function
loglikelihood <- function(par, path_l, path_n, path_l_cscore, path_l_oltv, path_l_dti, path_l_origamt, path_l_lastrt, tobs, cobs){
  f1 <- par[grepl("f1",names(par))]
  Zc <- par[grepl("Zc",names(par))]
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
  
  f_[1] <- f1
  
  #compute likelihood and other elements at every t
  for(i in 1:length(tobs)){
    sum_cov <- Beta[1]*path_l_cscore[DateQtr == tobs[i]]$CSCORE_B + Beta[2]*path_l_oltv[DateQtr == tobs[i]]$OLTV + 
      Beta[3]*path_l_dti[DateQtr == tobs[i]]$DTI + Beta[4]*path_l_origamt[DateQtr == tobs[i]]$ORIG_AMT +
      Beta[5]*path_l_lastrt[DateQtr == tobs[i]]$LAST_RT
    
    
    p_ <- 1/(1 + exp(Zc*f_[i] - sum_cov))
    
    score_l_[i] <- sum(p_)*Zc - cobs[i]*path_l[i]*Zc 
    score_n_[i] <- Zm%*%Siginv%*%(path_n[i,] - Zm*f_[i])
    
    score_[i] <-  score_l_[i] + score_n_[i]
    
    loglike_l <- -cobs[i]*path_l[i]*Zc*f_[i] - sum(log(1 + exp(sum_cov - Zc*f_[i])))
    loglike_n <- -0.5*ncol(path_n)*log(2*pi) - 0.5*log(det(Sig)) - 0.5 * t((path_n[i,] - Zm*f_[i])) %*% Siginv %*% (path_n[i,] - Zm*f_[i])  
    
    loglike[i] <- loglike_n + loglike_l
    
    f_[i+1] <- w + A*score_[i] + B*f_[i]
  }
  loglike <- sum(loglike)
  return(-loglike)
}

## Common Elements
#Logit Variables
path_l_cscore <- subset(DM_1C, select = c("CSCORE_B","DateQtr"))
path_l_oltv <- subset(DM_1C, select = c("OLTV","DateQtr"))
path_l_dti <- subset(DM_1C, select = c("DTI","DateQtr"))
path_l_origamt <- subset(DM_1C, select = c("ORIG_AMT","DateQtr"))
path_l_lastrt <- subset(DM_1C, select = c("LAST_RT","DateQtr"))


### 2003 - 2007
##Calibration
tobs <- unique(DM_1C[DateQtr >= 2003 & DateQtr < 2008]$DateQtr)
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

par_03_07_GAS <- par

#par <- c(par_03_07_GAS[1],  par_03_07_GAS[2], par_03_07_GAS[3], par_03_07_GAS[4], par_03_07_GAS[5], abs(par_03_07_GAS[6]), abs(par_03_07_GAS[7]), abs(par_03_07_GAS[8]), abs(par_03_07_GAS[9]),
#         par_03_07_GAS[10], par_03_07_GAS[11], par_03_07_GAS[12], par_03_07_GAS[13], par_03_07_GAS[14])


#Comparison of fitted default rate with observed default rate:
Beta_ <- c(par_03_07_GAS["Beta1"], par_03_07_GAS["Beta2"], par_03_07_GAS["Beta3"], par_03_07_GAS["Beta4"], par_03_07_GAS["Beta5"])
p_average <- matrix(data = NA, nrow = length(tobs), ncol = 1)
for(i in 1:length(tobs)){
  #Dynamic probability for logit component
  sum_cov <- Beta_[1]*path_l_cscore[DateQtr == tobs[i]]$CSCORE_B + Beta_[2]*path_l_oltv[DateQtr == tobs[i]]$OLTV + 
    Beta_[3]*path_l_dti[DateQtr == tobs[i]]$DTI + Beta_[4]*path_l_origamt[DateQtr == tobs[i]]$ORIG_AMT +
    Beta_[5]*path_l_lastrt[DateQtr == tobs[i]]$LAST_RT
  
  p_average[i] <- sum(1/(1 + exp(par_03_07_GAS["Zc1"]*f_[i] - sum_cov)))/cobs[i]
}

plot(p_average, type = "l")
plot(path_l, type = "l")

#Unconditional mean of Factor
Ef <- fit$par["w1"]/(1 - fit$par["B1"])

plot(score_, type = "l")
plot(f_, type = "l")
f_03_07 <- f_
score_03_07 <- score_


##Forecast
tobs <- unique(DM_1C[DateQtr >= 2003 & DateQtr < 2010]$DateQtr)
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

#***Run the iterative procedure within the likelihood function using the parameters par_03_07_GAS

DM_1C[, Factor := NA_real_]
#Assign factor
for(i in 1:length(tobs)){
  DM_1C[DateQtr == tobs[i], Factor := f_[i]]
}

#***Use parameters with 1% sig.
DM_1C[DateQtr >= 2008 & DateQtr < 2010, PD_forecast_GAS := exp(par["Beta1"]*CSCORE_B + par["Beta5"]*LAST_RT - 
                                                                 par["Zc1"]*Factor)/(1+ exp(par["Beta1"]  + 
                                                                                              par["Beta2"]*CSCORE_B + par["Beta5"]*LAST_RT - par["Zc1"]*Factor))]

##Performance 
#AUC
auc_03_07_GAS_1 <- auc(DM_1C[DateQtr == 2008.00]$Default, DM_1C[DateQtr == 2008.00]$PD_forecast_GAS)
auc_03_07_GAS_2 <- auc(DM_1C[DateQtr == 2008.25]$Default, DM_1C[DateQtr == 2008.25]$PD_forecast_GAS)
auc_03_07_GAS_3 <- auc(DM_1C[DateQtr == 2008.50]$Default, DM_1C[DateQtr == 2008.50]$PD_forecast_GAS)
auc_03_07_GAS_4 <- auc(DM_1C[DateQtr == 2008.75]$Default, DM_1C[DateQtr == 2008.75]$PD_forecast_GAS)
auc_03_07_GAS_5 <- auc(DM_1C[DateQtr == 2009.00]$Default, DM_1C[DateQtr == 2009.00]$PD_forecast_GAS)
auc_03_07_GAS_6 <- auc(DM_1C[DateQtr == 2009.25]$Default, DM_1C[DateQtr == 2009.25]$PD_forecast_GAS)
auc_03_07_GAS_7 <- auc(DM_1C[DateQtr == 2009.50]$Default, DM_1C[DateQtr == 2009.50]$PD_forecast_GAS)
auc_03_07_GAS_8 <- auc(DM_1C[DateQtr == 2009.75]$Default, DM_1C[DateQtr == 2009.75]$PD_forecast_GAS)

#MAD
tobs <- unique(DM_1C[DateQtr >= 2008 & DateQtr < 2010]$DateQtr)
def_rate <- matrix(data = NA, nrow = length(tobs), ncol =1)
def_rate_estim <- matrix(data = NA, nrow = length(tobs), ncol = 1)

for(t in 1:length(tobs)){
  def_rate_estim[t] <- sum(DM_1C[DateQtr == tobs[t]]$PD_forecast_GAS)/nrow(DM_1C[DateQtr == tobs[t]])
  def_rate[t] <- sum(DM_1C[DateQtr == tobs[t]]$Default)/nrow(DM_1C[DateQtr == tobs[t]])
}

MAD_03_07_GAS <- sum(abs(def_rate - def_rate_estim))/length(tobs)



### 2003 - 2014
##Calibration
tobs <- unique(DM_1C[DateQtr >= 2003 & DateQtr < 2014]$DateQtr)
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

par_03_14_GAS <- par

#par <- c(par_03_14_GAS[1],  par_03_14_GAS[2], par_03_14_GAS[3], par_03_14_GAS[4], par_03_14_GAS[5], abs(par_03_14_GAS[6]), abs(par_03_14_GAS[7]), abs(par_03_14_GAS[8]), abs(par_03_14_GAS[9]),
#         par_03_14_GAS[10], par_03_14_GAS[11], par_03_14_GAS[12], par_03_14_GAS[13], par_03_14_GAS[14])

#Comparison of fitted default rate with observed default rate:
Beta_ <- c(par_03_14_GAS["Beta1"], par_03_14_GAS["Beta2"], par_03_14_GAS["Beta3"], par_03_14_GAS["Beta4"], par_03_14_GAS["Beta5"])
p_average <- matrix(data = NA, nrow = length(tobs), ncol = 1)
for(i in 1:length(tobs)){
  #Dynamic probability for logit component
  sum_cov <- Beta_[1]*path_l_cscore[DateQtr == tobs[i]]$CSCORE_B + Beta_[2]*path_l_oltv[DateQtr == tobs[i]]$OLTV + 
    Beta_[3]*path_l_dti[DateQtr == tobs[i]]$DTI + Beta_[4]*path_l_origamt[DateQtr == tobs[i]]$ORIG_AMT +
    Beta_[5]*path_l_lastrt[DateQtr == tobs[i]]$LAST_RT
  
  p_average[i] <- sum(1/(1 + exp(par_03_14_GAS["Zc1"]*f_[i] - sum_cov)))/cobs[i]
}

plot(p_average, type = "l")
plot(path_l, type = "l")

#Unconditional mean of Factor
Ef <- fit$par["w1"]/(1 - fit$par["B1"])

plot(score_, type = "l")
plot(f_, type = "l")
f_03_14 <- f_
score_03_14 <- score_


##Forecast
tobs <- unique(DM_1C[DateQtr >= 2003 & DateQtr < 2016]$DateQtr)
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

#***Run the iterative procedure within the likelihood function using the parameters par_03_07_GAS

DM_1C[, Factor := NA_real_]
#Assign factor
for(i in 1:length(tobs)){
  DM_1C[DateQtr == tobs[i], Factor := f_[i]]
}

#***Use parameters with 1% sig.
DM_1C[DateQtr >= 2014 & DateQtr < 2016, PD_forecast_GAS := exp(par["Beta1"]*CSCORE_B  + 
                                                                 # par["Beta2"]*OLTV + par["Beta3"]*DTI + par["Beta4"]*ORIG_AMT +
                                                                 par["Beta5"]*LAST_RT -
                                                                 par["Zc1"]*Factor)/(1+ exp(par["Beta1"]*CSCORE_B  + 
                                                                                              #par["Beta2"]*OLTV + par["Beta3"]*DTI + par["Beta4"]*ORIG_AMT +  
                                                                                              par["Beta5"]*LAST_RT - par["Zc1"]*Factor))]

##Performance 
#AUC
auc_03_14_GAS_1 <- auc(DM_1C[DateQtr == 2014.00]$Default, DM_1C[DateQtr == 2014.00]$PD_forecast_GAS)
auc_03_14_GAS_2 <- auc(DM_1C[DateQtr == 2014.25]$Default, DM_1C[DateQtr == 2014.25]$PD_forecast_GAS)
auc_03_14_GAS_3 <- auc(DM_1C[DateQtr == 2014.50]$Default, DM_1C[DateQtr == 2014.50]$PD_forecast_GAS)
auc_03_14_GAS_4 <- auc(DM_1C[DateQtr == 2014.75]$Default, DM_1C[DateQtr == 2014.75]$PD_forecast_GAS)
auc_03_14_GAS_5 <- auc(DM_1C[DateQtr == 2015.00]$Default, DM_1C[DateQtr == 2015.00]$PD_forecast_GAS)
auc_03_14_GAS_6 <- auc(DM_1C[DateQtr == 2015.25]$Default, DM_1C[DateQtr == 2015.25]$PD_forecast_GAS)
auc_03_14_GAS_7 <- auc(DM_1C[DateQtr == 2015.50]$Default, DM_1C[DateQtr == 2015.50]$PD_forecast_GAS)
auc_03_14_GAS_8 <- auc(DM_1C[DateQtr == 2015.75]$Default, DM_1C[DateQtr == 2015.75]$PD_forecast_GAS)


#MAD
tobs <- unique(DM_1C[DateQtr >= 2014 & DateQtr < 2016]$DateQtr)
def_rate <- matrix(data = NA, nrow = length(tobs), ncol =1)
def_rate_estim <- matrix(data = NA, nrow = length(tobs), ncol = 1)

for(t in 1:length(tobs)){
  def_rate_estim[t] <- sum(DM_1C[DateQtr == tobs[t]]$PD_forecast_GAS)/nrow(DM_1C[DateQtr == tobs[t]])
  def_rate[t] <- sum(DM_1C[DateQtr == tobs[t]]$Default)/nrow(DM_1C[DateQtr == tobs[t]])
}

MAD_03_14_GAS <- sum(abs(def_rate - def_rate_estim))/length(tobs)



### Robustness check for different Scaling matrix ###

#Generate Series (same parameters!)
loglikelihood_GAS2 <- function(par, path_l, path_n, path_l_cscore, path_l_oltv, path_l_dti, path_l_origamt, path_l_lastrt, tobs, cobs){
  f1 <- par[grepl("f1",names(par))]
  Zc <- par[grepl("Zc",names(par))]
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
  
  f_[1] <- f1
  
  #compute likelihood and other elements at every t
  for(i in 1:length(tobs)){
    sum_cov <- Beta[1]*path_l_cscore[DateQtr == tobs[i]]$CSCORE_B + Beta[2]*path_l_oltv[DateQtr == tobs[i]]$OLTV + 
      Beta[3]*path_l_dti[DateQtr == tobs[i]]$DTI + Beta[4]*path_l_origamt[DateQtr == tobs[i]]$ORIG_AMT +
      Beta[5]*path_l_lastrt[DateQtr == tobs[i]]$LAST_RT
    
    
    p_ <- 1/(1 + exp(Zc*f_[i] - sum_cov))
    
    score_l_[i] <- sum(p_)*Zc - cobs[i]*path_l[i]*Zc 
    score_n_[i] <- Zm%*%Siginv%*%(path_n[i,] - Zm*f_[i])
    
    #Fisher Information Matrix
    sum_comp <- 0 
    for(n in 1:length(p_)){
      sum_comp <- sum_comp + t(p_)%*%p_ - p_[n]^2
    }
    I_c <- Zc%*%t(Zc)*as.numeric(sum(p_) + sum_comp - sum(p_)^2)
    I_m <- t(Zm)%*%Siginv%*%Zm
    I <- I_c + I_m
    
    #Scaling matrix
    ev <<- eigen(I)
    L  <- ev$values[which(ev$values != 0)]
    V <- as.matrix(ev$vectors[,which(ev$values != 0)])
    if(length(L) == 0){
      S <- matrix(data = 0, nrow = ncol(f_), ncol = ncol(f_))
    }else if(length(L) == ncol(V)){
      S <- V %*% diag(1/L, ncol = length(L), nrow = length(L))^(1/2) %*% t(V)
    }
    
    score_[i] <- S%*%(score_l_[i] + score_n_[i])
    
    loglike_l <- -cobs[i]*path_l[i]*Zc*f_[i] - sum(log(1 + exp(sum_cov - Zc*f_[i])))
    loglike_n <- -0.5*ncol(path_n)*log(2*pi) - 0.5*log(det(Sig)) - 0.5 * t((path_n[i,] - Zm*f_[i])) %*% Siginv %*% (path_n[i,] - Zm*f_[i])  
    
    loglike[i] <- loglike_n + loglike_l
    
    f_[i+1] <- w +  A*score_[i] + B*f_[i]
  }
  loglike <- sum(loglike)
  return(-loglike)
}


### 2003 - 2007
##Calibration
tobs <- unique(DM_1C[DateQtr >= 2003 & DateQtr < 2008]$DateQtr)
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

#Transform parameters to their restricted counter part
par <- c(par_03_07_GAS[1],  par_03_07_GAS[2], par_03_07_GAS[3], par_03_07_GAS[4], par_03_07_GAS[5], abs(par_03_07_GAS[6]), abs(par_03_07_GAS[7]), abs(par_03_07_GAS[8]), abs(par_03_07_GAS[9]),
         par_03_07_GAS[10], par_03_07_GAS[11], par_03_07_GAS[12], par_03_07_GAS[13], par_03_07_GAS[14])

#Comparison of fitted default rate with observed default rate:
Beta_ <- c(par_03_07_GAS["Beta1"], par_03_07_GAS["Beta2"],par_03_07_GAS["Beta3"], par_03_07_GAS["Beta4"], par_03_07_GAS["Beta5"])
p_average <- matrix(data = NA, nrow = length(tobs), ncol = 1)
for(i in 1:length(tobs)){
  #Dynamic probability for logit component
  sum_cov <- Beta_[1]*path_l_cscore[DateQtr == tobs[i]]$CSCORE_B + Beta_[2]*path_l_oltv[DateQtr == tobs[i]]$OLTV + 
    Beta_[3]*path_l_dti[DateQtr == tobs[i]]$DTI + Beta_[4]*path_l_origamt[DateQtr == tobs[i]]$ORIG_AMT +
    Beta_[5]*path_l_lastrt[DateQtr == tobs[i]]$LAST_RT
  
  p_average[i] <- sum(1/(1 + exp(par_03_07_GAS["Zc1"]*f_[i] - sum_cov)))/cobs[i]
}


plot(p_average, type = "l")
plot(path_l, type = "l")

#Unconditional mean of Factor
Ef <- fit$par["w1"]/(1 - fit$par["B1"])

plot(score_, type = "l")
plot(f_, type = "l")
f_03_07_GAS2 <- f_
score_03_07_GAS2 <- score_

##Forecast
tobs <- unique(DM_1C[DateQtr >= 2003 & DateQtr < 2010]$DateQtr)
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

#***Run the iterative procedure within the likelihood function using the parameters par_03_07_GAS

DM_1C[, Factor := NA_real_]
#Assign factor
for(i in 1:length(tobs)){
  DM_1C[DateQtr == tobs[i], Factor := f_[i]]
}

#***Use parameters with 1% sig.
DM_1C[DateQtr >= 2008 & DateQtr < 2010, PD_forecast_GAS := NA]
DM_1C[DateQtr >= 2008 & DateQtr < 2010, PD_forecast_GAS := exp(par["Beta1"]*CSCORE_B  + 
                                                                 # par["Beta2"]*OLTV + par["Beta3"]*DTI + par["Beta4"]*ORIG_AMT +
                                                                 par["Beta5"]*LAST_RT -
                                                                 par["Zc1"]*Factor)/(1+ exp(par["Beta1"]*CSCORE_B  + 
                                                                                              #par["Beta2"]*OLTV + par["Beta3"]*DTI + par["Beta4"]*ORIG_AMT +  
                                                                                              par["Beta5"]*LAST_RT - par["Zc1"]*Factor))]


auc_03_07_GAS2_1 <- auc(DM_1C[DateQtr == 2008.00]$Default, DM_1C[DateQtr == 2008.00]$PD_forecast_GAS)
auc_03_07_GAS2_2 <- auc(DM_1C[DateQtr == 2008.25]$Default, DM_1C[DateQtr == 2008.25]$PD_forecast_GAS)
auc_03_07_GAS2_3 <- auc(DM_1C[DateQtr == 2008.50]$Default, DM_1C[DateQtr == 2008.50]$PD_forecast_GAS)
auc_03_07_GAS2_4 <- auc(DM_1C[DateQtr == 2008.75]$Default, DM_1C[DateQtr == 2008.75]$PD_forecast_GAS)
auc_03_07_GAS2_5 <- auc(DM_1C[DateQtr == 2009.00]$Default, DM_1C[DateQtr == 2009.00]$PD_forecast_GAS)
auc_03_07_GAS2_6 <- auc(DM_1C[DateQtr == 2009.25]$Default, DM_1C[DateQtr == 2009.25]$PD_forecast_GAS)
auc_03_07_GAS2_7 <- auc(DM_1C[DateQtr == 2009.50]$Default, DM_1C[DateQtr == 2009.50]$PD_forecast_GAS)
auc_03_07_GAS2_8 <- auc(DM_1C[DateQtr == 2009.75]$Default, DM_1C[DateQtr == 2009.75]$PD_forecast_GAS)

#MAD
tobs <- unique(DM_1C[DateQtr >= 2008 & DateQtr < 2010]$DateQtr)
def_rate <- matrix(data = NA, nrow = length(tobs), ncol =1)
def_rate_estim <- matrix(data = NA, nrow = length(tobs), ncol = 1)

for(t in 1:length(tobs)){
  def_rate_estim[t] <- sum(DM_1C[DateQtr == tobs[t]]$PD_forecast_GAS)/nrow(DM_1C[DateQtr == tobs[t]])
  def_rate[t] <- sum(DM_1C[DateQtr == tobs[t]]$Default)/nrow(DM_1C[DateQtr == tobs[t]])
}

MAD_03_07_GAS2 <- sum(abs(def_rate - def_rate_estim))/length(tobs)


### 2003 - 2014
##Calibration
tobs <- unique(DM_1C[DateQtr >= 2003 & DateQtr < 2014]$DateQtr)
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

#Transform parameters to their restricted counter part
par <- c(par_03_14_GAS[1],  par_03_14_GAS[2], par_03_14_GAS[3], par_03_14_GAS[4], par_03_14_GAS[5], abs(par_03_14_GAS[6]), abs(par_03_14_GAS[7]), abs(par_03_14_GAS[8]), abs(par_03_14_GAS[9]),
         par_03_14_GAS[10], par_03_14_GAS[11], par_03_14_GAS[12], par_03_14_GAS[13], par_03_14_GAS[14])

#Comparison of fitted default rate with observed default rate:
Beta_ <- c(par_03_14_GAS["Beta1"], par_03_14_GAS["Beta2"],par_03_14_GAS["Beta3"], par_03_14_GAS["Beta4"], par_03_14_GAS["Beta5"])
p_average <- matrix(data = NA, nrow = length(tobs), ncol = 1)
for(i in 1:length(tobs)){
  #Dynamic probability for logit component
  sum_cov <- Beta_[1]*path_l_cscore[DateQtr == tobs[i]]$CSCORE_B + Beta_[2]*path_l_oltv[DateQtr == tobs[i]]$OLTV + 
    Beta_[3]*path_l_dti[DateQtr == tobs[i]]$DTI + Beta_[4]*path_l_origamt[DateQtr == tobs[i]]$ORIG_AMT +
    Beta_[5]*path_l_lastrt[DateQtr == tobs[i]]$LAST_RT
  
  p_average[i] <- sum(1/(1 + exp(par_03_14_GAS["Zc1"]*f_[i] - sum_cov)))/cobs[i]
}


plot(p_average, type = "l")
plot(path_l, type = "l")

#Unconditional mean of Factor
Ef <- fit$par["w1"]/(1 - fit$par["B1"])

plot(score_, type = "l")
plot(f_, type = "l")
f_03_07_GAS2 <- f_
score_03_07_GAS2 <- score_

##Forecast
tobs <- unique(DM_1C[DateQtr >= 2003 & DateQtr < 2016]$DateQtr)
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

#***Run the iterative procedure within the likelihood function using the parameters par_03_14_GAS

DM_1C[, Factor := NA_real_]
#Assign factor
for(i in 1:length(tobs)){
  DM_1C[DateQtr == tobs[i], Factor := f_[i]]
}

#***Use parameters with 1% sig.
DM_1C[DateQtr >= 2014 & DateQtr < 2016, PD_forecast_GAS := NA]
DM_1C[DateQtr >= 2014 & DateQtr < 2016, PD_forecast_GAS := exp(par["Beta1"]*CSCORE_B  + 
                                                                 # par["Beta2"]*OLTV + par["Beta3"]*DTI + par["Beta4"]*ORIG_AMT +
                                                                 par["Beta5"]*LAST_RT -
                                                                 par["Zc1"]*Factor)/(1+ exp(par["Beta1"]*CSCORE_B  + 
                                                                                              #par["Beta2"]*OLTV + par["Beta3"]*DTI + par["Beta4"]*ORIG_AMT +  
                                                                                              par["Beta5"]*LAST_RT - par["Zc1"]*Factor))]

##Performance 
#AUC
auc_03_14_GAS2_1 <- auc(DM_1C[DateQtr == 2014.00]$Default, DM_1C[DateQtr == 2014.00]$PD_forecast_GAS)
auc_03_14_GAS2_2 <- auc(DM_1C[DateQtr == 2014.25]$Default, DM_1C[DateQtr == 2014.25]$PD_forecast_GAS)
auc_03_14_GAS2_3 <- auc(DM_1C[DateQtr == 2014.50]$Default, DM_1C[DateQtr == 2014.50]$PD_forecast_GAS)
auc_03_14_GAS2_4 <- auc(DM_1C[DateQtr == 2014.75]$Default, DM_1C[DateQtr == 2014.75]$PD_forecast_GAS)
auc_03_14_GAS2_5 <- auc(DM_1C[DateQtr == 2015.00]$Default, DM_1C[DateQtr == 2015.00]$PD_forecast_GAS)
auc_03_14_GAS2_6 <- auc(DM_1C[DateQtr == 2015.25]$Default, DM_1C[DateQtr == 2015.25]$PD_forecast_GAS)
auc_03_14_GAS2_7 <- auc(DM_1C[DateQtr == 2015.50]$Default, DM_1C[DateQtr == 2015.50]$PD_forecast_GAS)
auc_03_14_GAS2_8 <- auc(DM_1C[DateQtr == 2015.75]$Default, DM_1C[DateQtr == 2015.75]$PD_forecast_GAS)


#MAD
tobs <- unique(DM_1C[DateQtr >= 2014 & DateQtr < 2016]$DateQtr)
def_rate <- matrix(data = NA, nrow = length(tobs), ncol =1)
def_rate_estim <- matrix(data = NA, nrow = length(tobs), ncol = 1)

for(t in 1:length(tobs)){
  def_rate_estim[t] <- sum(DM_1C[DateQtr == tobs[t]]$PD_forecast_GAS)/nrow(DM_1C[DateQtr == tobs[t]])
  def_rate[t] <- sum(DM_1C[DateQtr == tobs[t]]$Default)/nrow(DM_1C[DateQtr == tobs[t]])
}

MAD_03_14_GAS2 <- sum(abs(def_rate - def_rate_estim))/length(tobs)


