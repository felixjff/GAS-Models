#### Purpose: Determine the effect of the scaling matrix on the generated factors and success series. 

### Identity Matrix ###

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
loglikelihood_I <- function(par, path_l, path_n, tobs, cobs){
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
  scaled_score_ <- matrix(data = NA, nrow = length(tobs), ncol = 2)
  f_ <- matrix(data = NA, nrow = length(tobs)+1, ncol = 2)
  loglike <- matrix(data = NA, nrow = length(tobs), ncol = 1)
  #Logit
  p_ <- matrix(data = NA, nrow = length(tobs), ncol = 1)
  score_l_ <- matrix(data = NA, nrow = length(tobs), ncol = 2)
  #Normal
  score_n_ <- matrix(data = NA, nrow = length(tobs), ncol = 2)
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
    
    #Identity Scaling matrix
    S <- diag(1, ncol(f_), ncol(f_))
    
    #Scaled score
    scaled_score_[i,] <- S%*%score_[i,]
    
    #Log-likelihood
    loglike_l <- -cobs[i]*path_l[i]%*%Zc%*%f_[i,] - cobs[i]*log(1 + exp(-Zc%*%f_[i,]))
    loglike_n <- -0.5*ncol(path_n)*log(2*pi) - 0.5*log(det(Sig)) - 0.5 * t((path_n[i] - Zm%*%f_[i,])) %*% Siginv %*% (path_n[i] - Zm%*%f_[i,])  
    
    loglike[i] <- loglike_n + loglike_l
    
    f_[i+1,] <- w +  A%*%scaled_score_[i,] + B%*%f_[i,]
  }
  #compute log-likelihood
  loglike <- sum(loglike)
  return(-loglike)
}

#Initialization
parameters <- c(Zm1 = 0.5, Zc1 = 0.001, Zc2 = 0.001, B1 = 0.5, B2 = 0.5, A1 = 0.5, A2 = 0.5, Sig1 = 0.7, Sig2 = 0.1, f1m = 0, 
                f1c = 0 , w1 = 0, w2 = 0)

#Estimation
fit <- optim(par = parameters, fn = loglikelihood_I, method = "BFGS" , path_l = path_l, 
             path_n = path_n, cobs = cobs, tobs = tobs, control=list(trace = 1, REPORT=10, maxit = 300), hessian = TRUE)

#Transform parameters to their restricted counter part
par <- c(fit$par[1],  fit$par[2], fit$par[3], fit$par[4], fit$par[5], abs(fit$par[6]), abs(fit$par[7]), abs(fit$par[8]), abs(fit$par[9]),
         fit$par[10], fit$par[11], fit$par[12], fit$par[13])


#Comparison of fitted default rate with observed default rate:
p_I <- plot(p_, type = "l")
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

s_I1 <- plot(score_[,1], type = "l")
s_I2 <- plot(score_[,2], type = "l")
f_I1 <- plot(f_[,1], type = "l")
f_I2 <- plot(f_[,2], type = "l")


#### Squared Root Inverse of Eigendecomposition of Fisher Inormation Matrix ####

#Likelihood function
loglikelihood_Creal <- function(par, path_l, path_n, tobs, cobs){
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
  scaled_score_ <- matrix(data = NA, nrow = length(tobs), ncol = 2)
  f_ <- matrix(data = NA, nrow = length(tobs)+1, ncol = 2)
  loglike <- matrix(data = NA, nrow = length(tobs), ncol = 1)
  #Logit
  p_ <- matrix(data = NA, nrow = length(tobs), ncol = 1)
  score_l_ <- matrix(data = NA, nrow = length(tobs), ncol = 2)
  #Normal
  score_n_ <- matrix(data = NA, nrow = length(tobs), ncol = 2)
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
    
    #Fisher Information Matrix
    I_c <- cobs[i]*p_[i]*(1 - p_[i])*Zm%*%t(Zm)
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
    
    #Scaled score
    scaled_score_[i,] <- S%*%score_[i,]
    
    #Log-likelihood
    loglike_l <- -cobs[i]*path_l[i]%*%Zc%*%f_[i,] - cobs[i]*log(1 + exp(-Zc%*%f_[i,]))
    loglike_n <- -0.5*ncol(path_n)*log(2*pi) - 0.5*log(det(Sig)) - 0.5 * t((path_n[i] - Zm%*%f_[i,])) %*% Siginv %*% (path_n[i] - Zm%*%f_[i,])  
    
    loglike[i] <- loglike_n + loglike_l
    
    f_[i+1,] <- w + A%*%scaled_score_[i,] + B%*%f_[i,]
  }
  #compute log-likelihood
  loglike <- sum(loglike)
  return(-loglike)
}

#Initialization
parameters <- c(Zm1 = 0.5, Zc1 = 0.1, Zc2 = 0.1, B1 = 0.5, B2 = 0.5, A1 = 0.05, A2 = 0.05,  Sig1 = 0.4, Sig2 = 0.4, f1m = 0, 
                f1c = 0, w1 = 1, w2 = 1)

#Estimation
fit <- optim(par = parameters, fn = loglikelihood_Creal, method = "BFGS" , path_l = path_l, 
             path_n = path_n, cobs = cobs, tobs = tobs, control=list(trace = 1, REPORT=10, maxit = 300), hessian = TRUE)

#Transform parameters to their restricted counter part
  par <- c(fit$par[1],  fit$par[2], fit$par[3], fit$par[4], fit$par[5], abs(fit$par[6]), abs(fit$par[7]), abs(fit$par[8]), abs(fit$par[9]),
         fit$par[10], fit$par[11], fit$par[12], fit$par[13])


#Comparison of fitted default rate with observed default rate:
p_C <- plot(p_, type = "l")
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

s_C1 <- plot(score_[,1], type = "l")
s_C2 <- plot(score_[,2], type = "l")
f_C1 <- plot(f_[,1], type = "l")
f_C2 <- plot(f_[,2], type = "l")


#Conclusion:
#1. In this emprical exercise, the scaling matrix proposed by Creal et al. has been successfully introduced to scale the score(s).
#   The results on both sides are very clear, namely the scaling matrix does not have a significant effect on the estimation procedure. 
#2. In both scenarios, the maximization reached the same maximum and provided the same results for the generated rate series.
#3. The estimated factors are very different from each other. In both cases, one of the factors took the same form and the other ones
#   have different form. Furthermore, the starting value of one of the factors seems to be unrealistic in both cases is unrealistic. 
#   This is a clear hint towards overparametization and misspecification of the DGP. It seems that the second frailty factor is
#   unnecessary. 

