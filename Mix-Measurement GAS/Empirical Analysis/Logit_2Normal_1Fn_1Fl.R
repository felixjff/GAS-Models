#### Purpose: Calibrate the logit-normal GAS model with a pool of loans, 2 macro variables, one macro factor and one frailty factor

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
    p_[i] <- 1/(1 + exp(Zc%*%f_[i,]))
    
    #Score
    score_l_[i,] <- cobs[i]*p_[i]%*%Zc - cobs[i]*path_l[i]%*%Zc 
    score_n_[i,] <- Zm%*%Siginv%*%(path_n[i] - Zm%*%f_[i,])
    
    score_[i,] <-  score_l_[i,] + score_n_[i,]
    
    #Log-likelihood
    loglike_l <- -cobs[i]*path_l[i]%*%Zc%*%f_[i,] - cobs[i]*log(1 + exp(-Zc%*%f_[i,]))
    loglike_n <- -0.5*ncol(path_n)*log(2*pi) - 0.5*log(det(Sig)) - 0.5 * t((path_n[i] - Zm%*%f_[i,])) %*% Siginv %*% (path_n[i] - Zm%*%f_[i,])  
    
    loglike[i] <- loglike_n + loglike_l
    
    f_[i+1,] <- w +  A%*%score_[i,] + B%*%f_[i,]
  }
  #compute log-likelihood
  loglike <- sum(loglike)
  return(-loglike)
}

#Initialization
parameters <- c(Zm1 = 0.5, Zc1 = 0.001, Zc2 = 0.001, B1 = 0.5, B2 = 0.5, A1 = 0.5, A2 = 0.5, Sig1 = 0.7, Sig2 = 0.1, f1m = 0, 
                f1c = 0 , w1 = 0, w2 = 0)

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
Ef <- c(fit$par["w1"], fit$par["w2"])%*%solve((1 - diag(c(fit$par["B1"], fit$par["B2"]) , 2, 2)))

plot(score_[,1], type = "l")
plot(f_[,2], type = "l")

#Conclusions:
#1. The calibration results are extremely sensitive to the parameter constraints. 
#1.1 The parameter contrains on B to be within the unit interval generates a terrible p_ series. Although, the factors
#    and the scores are centered at their theoretical values.
#1.2 By only imposing the restriction on A, the generated p_ is very similar to the observed default rate. 
#1.3 By not imposing restrictions nor on A or B, the results are also terrible. 
#1.4 By setting the constraint w = 0 on the constant of the GAS process, the generated default rates seem to improve. 
#    However, the factors are not centred on their unconditional mean. The score is centred around zero.
#    The factors look more stable by using this constraint. 
#1.5 By imposing all three restrictions, the results are again terrible.
#1.6 The contraint on the diagonal of the factor loading matrix Zm has also shown to be essential for the maximisation to 
#    generate sensitive results. 