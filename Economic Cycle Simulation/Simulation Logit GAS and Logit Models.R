#Purpose: Simulation Study for Thesis. Based on binomial for every point in titme. 
#Author: Felix Farias

library(ggplot2)
library(stats4)
library(stats)
library(data.table)
#Simulate from normal distributions with defined mean and variance

####                              ####
#    Complete economic cycle      1  # 
####                              ####

sigma_nc <- 0.15
mu_nc <- log(7)
sigma_c <- 0.55
mu_c <- log(0.142)
sigma_r <- 0.55
mu_r <- log(0.777)

y_m_logit <- matrix(data = NA, nrow = 120, ncol = 10000)
sr_com_logit <- matrix(data = NA, nrow = 120, ncol = 10000) 
p_logit  <- matrix(data = NA, nrow = 120, ncol = 10000) 
tobs <- 120
cobs <- 1000
beta <- 1
tau <- 0.008

#1.common + 1.crisis + 1.recovery + 2.common + 2.crisis + 2.recovery
for(i in 1:10000){
  y_m_logit[1:36,i] <- rnorm(n = 36, mean = mu_nc, sd = sigma_nc)
  y_m_logit[37:48, i] <- rnorm(n = 6, mean = mu_c, sd = sigma_c)
  y_m_logit[49:56, i] <- rnorm(n = 8, mean = mu_r, sd = sigma_r)
  y_m_logit[57:92, i] <- rnorm(n = 36, mean = mu_nc, sd = sigma_nc)
  y_m_logit[93:104, i] <- rnorm(n = 12, mean = mu_c, sd = sigma_c)
  y_m_logit[105:120, i] <- rnorm(n = 16, mean = mu_r, sd = sigma_r)
  
  for(t in 1:tobs){
    #prob of default
    p_logit <- tau/(1+ exp(-beta*as.numeric(y_m_logit[t,i])))
    
    #Instead of storing all observations, only store the generated sucess rate.
    y_l <- rbinom(n = cobs, size = 1, prob = p_logit)
    sr_com_logit[t,i] <- sum(y_l)/cobs
  }
  print(i)
}

y_m_logit <- as.data.table(y_m_logit)
y_m_logit[, States := NA_character_]
y_m_logit[1:35, States := "Stable Growth"]
y_m_logit[36:47, States := "Crisis"]
y_m_logit[48:57, States := "Recovery"]
y_m_logit[58:91, States := "Stable Growth"]
y_m_logit[92:103, States := "Crisis"]
y_m_logit[104:120, States := "Recovery"]

#1.common + 1.crisis + 1.recovery + 2.common + 2.crisis + 2.recovery

par(mfrow = c(1,2))

gdp_p <- ggplot(simulations, aes(x = 1:120, y = V1, group = 1, colour = States)) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))  +
geom_line() + 
xlab("Observation") + 
ylab("Risk Factor Values")

pd_p <- ggplot(simulations, aes(x = 1:120, y = V501, group = 1, colour = States)) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))  +
  geom_line() + 
  xlab("Observation") + 
  ylab("Default Rate Values")

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

multiplot(gdp_p, pd_p)


###                              ####
#  Estimation with Given size of Cross-section      # 
####                             ####

#Simulation in 3 different setting varying the amount of obs. in the cross-section:

  Loglikelihood <- function(beta, path_n, path_l, tobs){
    cobs <- 1000
    tau <- 0.008
    loglike <- matrix(data = NA, nrow = tobs, ncol = 1)
    for(j in 1:tobs){
      p_ <- tau/(1+exp(-beta*as.numeric(path_n[j])))
      loglike[j] <- cobs*path_l[j]*beta*as.numeric(path_n[j]) - cobs*log(1 + exp(beta*as.numeric(path_n[j])))
    }
    
    
    logl <- sum(loglike)
    return(-logl)
  }

#Case 1: 1000 cross section and COMPLETE CYCLE 
beta <- 1
pe_com_logit <- matrix(data = NA, nrow = 10000, ncol = 1)
y_m_logit <- as.matrix(y_m_logit)

for(i in 1:10000){
  fit_model <- optim(par = beta, fn = Loglikelihood, method = "BFGS", 
                     path_n = y_m_logit[,i], path_l = sr_com_logit[,i], tobs = tobs)
  pe_com_logit[i] <- fit_model$par
  print(paste(i, pe_com_logit[i]))
}

bias_com_logit <- pe_com_logit - 1



#####
#####  2C, 2R
#####

#2C, 2R

y_m_2c2r_logit <- as.matrix(y_m_logit[States == "Crisis" | States == "Recovery"])
tobs <- nrow(y_m_2c2r_logit)
sr_2c2r_logit <- matrix(data = NA, nrow = tobs, ncol = 10000)

for(i in 1:10000){
  for(t in 1:tobs){
    #prob of default
    p_logit <- tau/(1+ exp(-beta*as.numeric(y_m_2c2r_logit[t,i])))
    
    #Instead of storing all observations, only store the generated sucess rate.
    y_l <- rbinom(n = cobs, size = 1, prob = p_logit)
    sr_2c2r_logit[t,i] <- sum(y_l)/cobs
  }
  print(i)
}

pe_2c2r_logit <- matrix(data = NA, nrow = 10000, ncol = 1)
y_m_2c2r_logit <- as.matrix(y_m_2c2r_logit)

for(i in 1:10000){
  fit_model <- optim(par = beta, fn = Loglikelihood, method = "BFGS", 
                     path_n = y_m_2c2r_logit[,i], path_l = sr_2c2r_logit[,i], tobs = tobs)
  pe_2c2r_logit[i] <- fit_model$par
  print(paste(i, pe_2c2r_logit[i]))
}

bias_2c2r_logit <- pe_2c2r_logit - 1

#####
#####  2C, 2SG
#####

y_m_2c2sg_logit <- as.matrix(y_m_logit[States == "Crisis" | States == "Stable Growth"])
tobs <- nrow(y_m_2c2sg_logit)
sr_2c2sg_logit <- matrix(data = NA, nrow = tobs, ncol = 10000)

for(i in 1:10000){
  for(t in 1:tobs){
    #prob of default
    p_logit <- tau/(1+ exp(-beta*as.numeric(y_m_2c2sg_logit[t,i])))
    
    #Instead of storing all observations, only store the generated sucess rate.
    y_l <- rbinom(n = cobs, size = 1, prob = p_logit)
    sr_2c2sg_logit[t,i] <- sum(y_l)/cobs
  }
  print(i)
}

pe_2c2sg_logit <- matrix(data = NA, nrow = 10000, ncol = 1)
y_m_2c2sg_logit <- as.matrix(y_m_2c2sg_logit)

for(i in 1:10000){
  fit_model <- optim(par = beta, fn = Loglikelihood, method = "BFGS", 
                     path_n = y_m_2c2sg_logit[,i], path_l = sr_2c2sg_logit[,i], tobs = tobs)
  pe_2c2sg_logit[i] <- fit_model$par
  print(paste(i, pe_2c2sg_logit[i]))
}

bias_2c2sg_logit <- pe_2c2sg_logit - 1

#####
#####  2R, 2SG
#####

y_m_2r2sg_logit <- as.matrix(y_m_logit[States == "Recovery" | States == "Stable Growth"])
tobs <- nrow(y_m_2r2sg_logit)
sr_2r2sg_logit <- matrix(data = NA, nrow = tobs, ncol = 10000)

for(i in 1:10000){
  for(t in 1:tobs){
    #prob of default
    p_logit <- tau/(1+ exp(-beta*as.numeric(y_m_2r2sg_logit[t,i])))
    
    #Instead of storing all observations, only store the generated sucess rate.
    y_l <- rbinom(n = cobs, size = 1, prob = p_logit)
    sr_2r2sg_logit[t,i] <- sum(y_l)/cobs
  }
  print(i)
}

pe_2r2sg_logit <- matrix(data = NA, nrow = 10000, ncol = 1)
y_m_2r2sg_logit <- as.matrix(y_m_2r2sg_logit)

for(i in 1:10000){
  fit_model <- optim(par = beta, fn = Loglikelihood, method = "BFGS", 
                     path_n = y_m_2r2sg_logit[,i], path_l = sr_2r2sg_logit[,i], tobs = tobs)
  pe_2r2sg_logit[i] <- fit_model$par
  print(paste(i, pe_2r2sg_logit[i]))
}

bias_2r2sg_logit <- pe_2r2sg_logit - 1

#####
#####  2C
#####

y_m_2c_logit <- as.matrix(y_m_logit[States == "Crisis"])
tobs <- nrow(y_m_2c_logit)
sr_2c_logit <- matrix(data = NA, nrow = tobs, ncol = 10000)

for(i in 1:10000){
  for(t in 1:tobs){
    #prob of default
    p_logit <- tau/(1+ exp(-beta*as.numeric(y_m_2c_logit[t,i])))
    
    #Instead of storing all observations, only store the generated sucess rate.
    y_l <- rbinom(n = cobs, size = 1, prob = p_logit)
    sr_2c_logit[t,i] <- sum(y_l)/cobs
  }
  print(i)
}

pe_2c_logit <- matrix(data = NA, nrow = 10000, ncol = 1)
y_m_2c_logit <- as.matrix(y_m_2c_logit)

for(i in 1:10000){
  fit_model <- optim(par = beta, fn = Loglikelihood, method = "BFGS", 
                     path_n = y_m_2c_logit[,i], path_l = sr_2c_logit[,i], tobs = tobs)
  pe_2c_logit[i] <- fit_model$par
  print(paste(i, pe_2c_logit[i]))
}

bias_2c_logit <- pe_2c_logit - 1

#####
#####  2R
#####

y_m_2r_logit <- as.matrix(y_m_logit[States == "Recovery"])
tobs <- nrow(y_m_2r_logit)
sr_2r_logit <- matrix(data = NA, nrow = tobs, ncol = 10000)

for(i in 1:10000){
  for(t in 1:tobs){
    #prob of default
    p_logit <- tau/(1+ exp(-beta*as.numeric(y_m_2r_logit[t,i])))
    
    #Instead of storing all observations, only store the generated sucess rate.
    y_l <- rbinom(n = cobs, size = 1, prob = p_logit)
    sr_2r_logit[t,i] <- sum(y_l)/cobs
  }
  print(i)
}

pe_2r_logit <- matrix(data = NA, nrow = 10000, ncol = 1)
y_m_2r_logit <- as.matrix(y_m_2r_logit)

for(i in 1:10000){
  fit_model <- optim(par = beta, fn = Loglikelihood, method = "BFGS", 
                     path_n = y_m_2r_logit[,i], path_l = sr_2r_logit[,i], tobs = tobs)
  pe_2r_logit[i] <- fit_model$par
  print(paste(i, pe_2r_logit[i]))
}

bias_2r_logit <- pe_2r_logit - 1


#####
#####  2SG
#####


y_m_2sg_logit <- as.matrix(y_m_logit[States == "Stable Growth"])
tobs <- nrow(y_m_2sg_logit)
sr_2sg_logit <- matrix(data = NA, nrow = tobs, ncol = 10000)

for(i in 1:10000){
  for(t in 1:tobs){
    #prob of default
    p_logit <- tau/(1+ exp(-beta*as.numeric(y_m_2sg_logit[t,i])))
    
    #Instead of storing all observations, only store the generated sucess rate.
    y_l <- rbinom(n = cobs, size = 1, prob = p_logit)
    sr_2sg_logit[t,i] <- sum(y_l)/cobs
  }
  print(i)
}

pe_2sg_logit <- matrix(data = NA, nrow = 10000, ncol = 1)
y_m_2sg_logit <- as.matrix(y_m_2sg_logit)

for(i in 1:10000){
  fit_model <- optim(par = beta, fn = Loglikelihood, method = "BFGS", 
                     path_n = y_m_2sg_logit[,i], path_l = sr_2sg_logit[,i], tobs = tobs)
  pe_2sg_logit[i] <- fit_model$par
  print(paste(i, pe_2sg_logit[i]))
}

bias_2sg_logit <- pe_2sg_logit - 1


#Tests for significant difference in bias -> Compare all to Complte sample
t.test(bias_com_logit, bias_2c2r_logit, alternative = "two.sided", conf.level = 0.95)
t.test(bias_com_logit, bias_2c2sg_logit, alternative = "two.sided", conf.level = 0.95)
t.test(bias_com_logit, bias_2r2sg_logit, alternative = "two.sided", conf.level = 0.95)
t.test(bias_com_logit, bias_2c_logit, alternative = "two.sided", conf.level = 0.95)
t.test(bias_com_logit, bias_2r_logit, alternative = "two.sided", conf.level = 0.95)
t.test(bias_com_logit, bias_2sg_logit, alternative = "two.sided", conf.level = 0.95)


#Quantifying the Bias
bias_2sg_logit <- as.data.table(bias_2sg_logit)
bias_2sg_logit[, SG := 1]
bias_2sg_logit[, C := 0]
bias_2sg_logit[, R := 0]
bias_2sg_logit[, TwoS := 0]
bias_2sg_logit[, Complete := 0]

bias_2r_logit <- as.data.table(bias_2r_logit)
bias_2r_logit[, SG := 0]
bias_2r_logit[, C := 0]
bias_2r_logit[, R := 1]
bias_2r_logit[, TwoS := 0]
bias_2r_logit[, Complete := 0]

bias_2c_logit <- as.data.table(bias_2c_logit)
bias_2c_logit[, SG := 0]
bias_2c_logit[, C := 1]
bias_2c_logit[, R := 0]
bias_2c_logit[, TwoS := 0]
bias_2c_logit[, Complete := 0]

bias_2r2sg_logit <- as.data.table(bias_2r2sg_logit)
bias_2r2sg_logit[, SG := 1]
bias_2r2sg_logit[, C := 0]
bias_2r2sg_logit[, R := 1]
bias_2r2sg_logit[, TwoS := 1]
bias_2r2sg_logit[, Complete := 0]

bias_2c2sg_logit <- as.data.table(bias_2c2sg_logit)
bias_2c2sg_logit[, SG := 1]
bias_2c2sg_logit[, C := 1]
bias_2c2sg_logit[, R := 0]
bias_2c2sg_logit[, TwoS := 1]
bias_2c2sg_logit[, Complete := 0]

bias_2c2r_logit <- as.data.table(bias_2c2r_logit)
bias_2c2r_logit[, SG := 0]
bias_2c2r_logit[, C := 1]
bias_2c2r_logit[, R := 1]
bias_2c2r_logit[, TwoS := 1]
bias_2c2r_logit[, Complete := 0]

bias_com_logit <- as.data.table(bias_com_logit)
bias_com_logit[, SG := 1]
bias_com_logit[, C := 1]
bias_com_logit[, R := 1]
bias_com_logit[, TwoS := 0]
bias_com_logit[, Complete := 0]

bias_all_logit <- rbind(bias_com_logit, bias_2r2sg_logit, bias_2c2r_logit, bias_2c2sg_logit, bias_2sg_logit, bias_2r_logit, bias_2c_logit)
quant_bias <- lm(V1 ~ SG + C + R, data = bias_all_logit)

#Previous to using this function save the data collected so far
summary.lm(quant_bias)





#############
#############   CREAL ET AL MODEL SIMULATION   ################
#############
  
  ###SIMULATION GAS MODEL###
  simulations <- 2000
  cobs <- 1000
  
  #GAS process
  score_com <- matrix(data = NA, nrow = tobs, ncol = simulations)
  f_com <- matrix(data = NA, nrow = tobs+1, ncol = simulations)
  #Logit
  p_com <- matrix(data = NA, nrow = tobs, ncol = simulations)
  gen_success_rate_com <- matrix(data = NA, nrow = tobs, ncol = simulations)
  y_m <- matrix(data = NA, nrow = tobs ,  ncol = simulations)
  
  #Normal - Economic Cycle
  sigma_nc <- 0.15
  mu_nc <- log(7)
  sigma_c <- 0.55
  mu_c <- log(0.142)
  sigma_r <- 0.55
  mu_r <- log(0.777)
  
  
  # y_m <- matrix(data = NA, nrow = tobs ,  ncol = simulations)
  # sig_nc <- 0.5
  # mu_nc <- log(12)
  # sig_c <- 1
  # mu_c <- log(0.00001)
  # sig_r <- 2
  # mu_r <- log(0.1)
  
  # Data Generating Process True Parameter Values
  true_param <- c(f1 = 0 , w = 2, alpha = 0.1, beta = 0.8, lambda = 0.001, delta = 0.3)
  f1 = 0 ; w = 0.08; alpha = 0.5; beta = 0.5; lambda = 0.001; delta = 0.3
  
  #Level 
  tau = 0.008
  
  
  #COMPLETE CYCLE
  tobs <- 120
  
  #simulate
  for(i in 1:simulations){
    
    y_m[1:36,i]    <- rnorm(n = 36, mean = mu_nc, sd = sigma_nc)
    y_m[37:48, i]  <- rnorm(n = 6, mean = mu_c, sd = sigma_c)
    y_m[49:56, i]  <- rnorm(n = 8, mean = mu_r, sd = sigma_r)
    y_m[57:92, i]  <- rnorm(n = 36, mean = mu_nc, sd = sigma_nc)
    y_m[93:104, i] <- rnorm(n = 12, mean = mu_c, sd = sigma_c)
    y_m[105:120, i]<- rnorm(n = 16, mean = mu_r, sd = sigma_r)
    
    #generate path
    f_com[1,i] <- true_param["f1"]
    for(t in 1:tobs){
      #Dynamic probability for logit component
      p_s <- true_param["w"] +  true_param["alpha"]*f_com[t,i] + true_param["beta"]*y_m[t,i]
      p_com[t,i] <- tau/(1+ exp(-p_s))
      
      #Instead of storing all observations, only store the generated sucess rate.
      y_l <- rbinom(n = cobs, size = 1, prob = p_com[t,i])
      gen_success_rate_com[t,i] <- sum(y_l)/cobs
      
      #scaling comp.
      #c <- 1/(1 - tau + exp( -true_param["w"] - true_param["alpha"]*f_com[t,i] - true_param["beta"]*y_m[t,i]))
      
      #Compute score
      #score_com[t,i] <- (cobs*gen_success_rate_com[t,i] - cobs*p_com[t,i]*c)*true_param["alpha"]
      score_com[t,i] <- (cobs*gen_success_rate_com[t,i] * (exp(-p_s)/(1-tau+exp(-p_s))) + 
                           ((1-tau)*exp(p_s))/(1+(1-tau)*exp(p_s)) - exp(p_s)/(1+exp(p_s)))*true_param["alpha"]
      
      
      f_com[t+1,i] <-  true_param["lambda"]*score_com[t,i] + true_param["delta"]*f_com[t,i]
    }
    print(i)
  }
  
  #For selection of cycle components
  y_m <- as.data.table(y_m)
  y_m[, States := NA_character_]
  y_m[1:35, States := "Stable Growth"]
  y_m[36:47, States := "Crisis"]
  y_m[48:57, States := "Recovery"]
  y_m[58:91, States := "Stable Growth"]
  y_m[92:103, States := "Crisis"]
  y_m[104:120, States := "Recovery"]

  
  
  #2C, 2R
  y_m_2c2r <- as.matrix(y_m[States == "Crisis" | States == "Recovery"])
  tobs <- nrow(y_m_2c2r)
  score_2c2r <- matrix(data = NA, nrow = tobs, ncol = simulations)
  f_2c2r <- matrix(data = NA, nrow = tobs+1, ncol = simulations)
  p_2c2r <- matrix(data = NA, nrow = tobs, ncol = simulations)
  gen_success_rate_2c2r <- matrix(data = NA, nrow = tobs, ncol = simulations)
  
for(i in 1:simulations){
  #generate path
  f_2c2r[1,i] <- true_param["f1"]
  for(t in 1:tobs){
    #Dynamic probability for logit component
    p_s <- true_param["w"] +  true_param["alpha"]*f_2c2r[t,i] + true_param["beta"]*as.numeric(y_m_2c2r[t,i])
    p_2c2r[t,i] <- tau/(1+ exp(-p_s))
    
    #Instead of storing all observations, only store the generated sucess rate.
    y_l <- rbinom(n = cobs, size = 1, prob = p_2c2r[t,i])
    gen_success_rate_2c2r[t,i] <- sum(y_l)/cobs
    
    #scaling comp.
    #c <- 1/(1 - tau + exp( -true_param["w"] - true_param["alpha"]*f_com[t,i] - true_param["beta"]*y_m[t,i]))
    
    #Compute score
    #score_com[t,i] <- (cobs*gen_success_rate_com[t,i] - cobs*p_com[t,i]*c)*true_param["alpha"]
    score_2c2r[t,i] <- (cobs*gen_success_rate_2c2r[t,i] * (exp(-p_s)/(1-tau+exp(-p_s))) + 
                         ((1-tau)*exp(p_s))/(1+(1-tau)*exp(p_s)) - exp(p_s)/(1+exp(p_s)))*true_param["alpha"]
    
    
    f_2c2r[t+1,i] <-  true_param["lambda"]*score_2c2r[t,i] + true_param["delta"]*f_2c2r[t,i]
  }
}

  
  #2C, 2SG
  y_m_2c2sg <- as.matrix(y_m[States == "Crisis" | States == "Stable Growth"])
  tobs <- nrow(y_m_2c2sg)
  score_2c2sg <- matrix(data = NA, nrow = tobs, ncol = simulations)
  f_2c2sg <- matrix(data = NA, nrow = tobs+1, ncol = simulations)
  p_2c2sg <- matrix(data = NA, nrow = tobs, ncol = simulations)
  gen_success_rate_2c2sg <- matrix(data = NA, nrow = tobs, ncol = simulations)
  
  for(i in 1:simulations){
    #generate path
    f_2c2sg[1,i] <- true_param["f1"]
    for(t in 1:tobs){
      #Dynamic probability for logit component
      p_s <- true_param["w"] +  true_param["alpha"]*f_2c2sg[t,i] + true_param["beta"]*as.numeric(y_m_2c2sg[t,i])
      p_2c2sg[t,i] <- tau/(1+ exp(-p_s))
      
      #Instead of storing all observations, only store the generated sucess rate.
      y_l <- rbinom(n = cobs, size = 1, prob = p_2c2sg[t,i])
      gen_success_rate_2c2sg[t,i] <- sum(y_l)/cobs
      
      #scaling comp.
      #c <- 1/(1 - tau + exp( -true_param["w"] - true_param["alpha"]*f_com[t,i] - true_param["beta"]*y_m[t,i]))
      
      #Compute score
      #score_com[t,i] <- (cobs*gen_success_rate_com[t,i] - cobs*p_com[t,i]*c)*true_param["alpha"]
      score_2c2sg[t,i] <- (cobs*gen_success_rate_2c2sg[t,i] * (exp(-p_s)/(1-tau+exp(-p_s))) + 
                            ((1-tau)*exp(p_s))/(1+(1-tau)*exp(p_s)) - exp(p_s)/(1+exp(p_s)))*true_param["alpha"]
      
      
      f_2c2sg[t+1,i] <-  true_param["lambda"]*score_2c2sg[t,i] + true_param["delta"]*f_2c2sg[t,i]
    }
  }
  
  

  #2R, 2SG
  y_m_2r2sg <- as.matrix(y_m[States == "Recovery" | States == "Stable Growth"])
  tobs <- nrow(y_m_2r2sg)
  score_2r2sg <- matrix(data = NA, nrow = tobs, ncol = simulations)
  f_2r2sg <- matrix(data = NA, nrow = tobs+1, ncol = simulations)
  p_2r2sg <- matrix(data = NA, nrow = tobs, ncol = simulations)
  gen_success_rate_2r2sg <- matrix(data = NA, nrow = tobs, ncol = simulations)
  
  for(i in 1:simulations){
    #generate path
    f_2r2sg[1,i] <- true_param["f1"]
    for(t in 1:tobs){
      p_s <- true_param["w"] +  true_param["alpha"]*f_2r2sg[t,i] + true_param["beta"]*as.numeric(y_m_2r2sg[t,i])
      p_2r2sg[t,i] <- tau/(1+ exp(-p_s))
      
      #Instead of storing all observations, only store the generated sucess rate.
      y_l <- rbinom(n = cobs, size = 1, prob = p_2r2sg[t,i])
      gen_success_rate_2r2sg[t,i] <- sum(y_l)/cobs
      
      #scaling comp.
      #c <- 1/(1 - tau + exp( -true_param["w"] - true_param["alpha"]*f_com[t,i] - true_param["beta"]*y_m[t,i]))
      
      #Compute score
      #score_com[t,i] <- (cobs*gen_success_rate_com[t,i] - cobs*p_com[t,i]*c)*true_param["alpha"]
      score_2r2sg[t,i] <- (cobs*gen_success_rate_2r2sg[t,i] * (exp(-p_s)/(1-tau+exp(-p_s))) + 
                             ((1-tau)*exp(p_s))/(1+(1-tau)*exp(p_s)) - exp(p_s)/(1+exp(p_s)))*true_param["alpha"]
      
      
      f_2r2sg[t+1,i] <-  true_param["lambda"]*score_2r2sg[t,i] + true_param["delta"]*f_2r2sg[t,i]
    }
  }
  
  
  #2C
  y_m_2c <- as.matrix(y_m[States == "Crisis"])
  tobs <- nrow(y_m_2c)
  score_2c <- matrix(data = NA, nrow = tobs, ncol = simulations)
  f_2c <- matrix(data = NA, nrow = tobs+1, ncol = simulations)
  p_2c <- matrix(data = NA, nrow = tobs, ncol = simulations)
  gen_success_rate_2c <- matrix(data = NA, nrow = tobs, ncol = simulations)
  
  for(i in 1:simulations){
    #generate path
    f_2c[1,i] <- true_param["f1"]
    for(t in 1:tobs){
      p_s <- true_param["w"] +  true_param["alpha"]*f_2c[t,i] + true_param["beta"]*as.numeric(y_m_2c[t,i])
      p_2c[t,i] <- tau/(1+ exp(-p_s))
      
      #Instead of storing all observations, only store the generated sucess rate.
      y_l <- rbinom(n = cobs, size = 1, prob = p_2c[t,i])
      gen_success_rate_2c[t,i] <- sum(y_l)/cobs
      
      #scaling comp.
      #c <- 1/(1 - tau + exp( -true_param["w"] - true_param["alpha"]*f_com[t,i] - true_param["beta"]*y_m[t,i]))
      
      #Compute score
      #score_com[t,i] <- (cobs*gen_success_rate_com[t,i] - cobs*p_com[t,i]*c)*true_param["alpha"]
      score_2c[t,i] <- (cobs*gen_success_rate_2c[t,i] * (exp(-p_s)/(1-tau+exp(-p_s))) + 
                             ((1-tau)*exp(p_s))/(1+(1-tau)*exp(p_s)) - exp(p_s)/(1+exp(p_s)))*true_param["alpha"]
      
      
      f_2c[t+1,i] <-  true_param["lambda"]*score_2c[t,i] + true_param["delta"]*f_2c[t,i]
    }
  }
  
  
  #2R
  y_m_2r <- as.matrix(y_m[States == "Recovery"])
  tobs <- nrow(y_m_2r)
  score_2r <- matrix(data = NA, nrow = tobs, ncol = simulations)
  f_2r <- matrix(data = NA, nrow = tobs+1, ncol = simulations)
  p_2r <- matrix(data = NA, nrow = tobs, ncol = simulations)
  gen_success_rate_2r <- matrix(data = NA, nrow = tobs, ncol = simulations)
  
  for(i in 1:simulations){
    #generate path
    f_2r[1,i] <- true_param["f1"]
    for(t in 1:tobs){
      p_s <- true_param["w"] +  true_param["alpha"]*f_2r[t,i] + true_param["beta"]*as.numeric(y_m_2r[t,i])
      p_2r[t,i] <- tau/(1+ exp(-p_s))
      
      #Instead of storing all observations, only store the generated sucess rate.
      y_l <- rbinom(n = cobs, size = 1, prob = p_2r[t,i])
      gen_success_rate_2r[t,i] <- sum(y_l)/cobs
      
      #scaling comp.
      #c <- 1/(1 - tau + exp( -true_param["w"] - true_param["alpha"]*f_com[t,i] - true_param["beta"]*y_m[t,i]))
      
      #Compute score
      #score_com[t,i] <- (cobs*gen_success_rate_com[t,i] - cobs*p_com[t,i]*c)*true_param["alpha"]
      score_2r[t,i] <- (cobs*gen_success_rate_2r[t,i] * (exp(-p_s)/(1-tau+exp(-p_s))) + 
                             ((1-tau)*exp(p_s))/(1+(1-tau)*exp(p_s)) - exp(p_s)/(1+exp(p_s)))*true_param["alpha"]
      
      
      f_2r[t+1,i] <-  true_param["lambda"]*score_2r[t,i] + true_param["delta"]*f_2r[t,i]
    }
  }
  
  
  #2SG
  y_m_2sg <- as.matrix(y_m[States == "Stable Growth"])
  tobs <- nrow(y_m_2sg)
  score_2sg <- matrix(data = NA, nrow = tobs, ncol = simulations)
  f_2sg <- matrix(data = NA, nrow = tobs+1, ncol = simulations)
  p_2sg <- matrix(data = NA, nrow = tobs, ncol = simulations)
  gen_success_rate_2sg <- matrix(data = NA, nrow = tobs, ncol = simulations)
  
  for(i in 1:simulations){
    #generate path
    f_2sg[1,i] <- true_param["f1"]
    for(t in 1:tobs){
      p_s <- true_param["w"] +  true_param["alpha"]*f_2sg[t,i] + true_param["beta"]*as.numeric(y_m_2sg[t,i])
      p_2sg[t,i] <- tau/(1+ exp(-p_s))
      
      #Instead of storing all observations, only store the generated sucess rate.
      y_l <- rbinom(n = cobs, size = 1, prob = p_2sg[t,i])
      gen_success_rate_2sg[t,i] <- sum(y_l)/cobs
      
      #scaling comp.
      #c <- 1/(1 - tau + exp( -true_param["w"] - true_param["alpha"]*f_com[t,i] - true_param["beta"]*y_m[t,i]))
      
      #Compute score
      #score_com[t,i] <- (cobs*gen_success_rate_com[t,i] - cobs*p_com[t,i]*c)*true_param["alpha"]
      score_2sg[t,i] <- (cobs*gen_success_rate_2sg[t,i] * (exp(-p_s)/(1-tau+exp(-p_s))) + 
                          ((1-tau)*exp(p_s))/(1+(1-tau)*exp(p_s)) - exp(p_s)/(1+exp(p_s)))*true_param["alpha"]
      
      
      f_2sg[t+1,i] <-  true_param["lambda"]*score_2sg[t,i] + true_param["delta"]*f_2sg[t,i]
    }
  }
  
  
  ###Log Likelihood Function###
  loglikelihood <- function(par, path_l, path_n, tobs, cobs){
    f1 <- par[grepl("f1",names(par))]
    w <- par[grepl("w",names(par))]
    lambda <- abs(par[grepl("lambda",names(par))])
    delta <- par[grepl("delta",names(par))]
    alpha <- par[grepl("alpha",names(par))]
    beta <- par[grepl("beta",names(par))]
    tau <- 0.008
    
    #Common
    score_ <- matrix(data = NA, nrow = tobs, ncol = 1)
    f_ <- matrix(data = NA, nrow = tobs+1, ncol = 1)
    loglike <- matrix(data = NA, nrow = tobs, ncol = 1)
    #Logit
    p_ <- matrix(data = NA, nrow = tobs, ncol = 1)
    
    f_[1] <- f1
    
    for(i in 1:tobs){
      #Estimated PD
      p_s <- w + alpha*f_[i] + beta*path_n[i]
      p_[i] <- tau/(1+ exp(-p_s))
      
      #scaling comp.
      #c <- 1/(1 - tau + exp( - w - alpha*f_[i] - beta*as.numeric(path_n[i])))
      
      
      #Score
      #score_[i] <- (cobs*path_l[i] - c*cobs*p_[i])*alpha
      score_[i] <- (cobs*path_l[i] * (exp(-p_s)/(1-tau+exp(-p_s))) + 
                           ((1-tau)*exp(p_s))/(1+(1-tau)*exp(p_s)) - exp(p_s)/(1+exp(p_s)))*alpha
    
      #Log-likelihood
      #loglike[i] <- cobs*path_l[i]*(w + alpha*f_[i] + beta*as.numeric(path_n[i])) - log(1 + p_[i])
      #loglike[i] <- log(1 + (1-tau)*exp(p_s)) - log(1+exp(p_s)) + cobs*path_l[i]*log(tau) -
      #              cobs*path_l[i]*log(1-tau+exp(-p_s))
      loglike[i] <- log(1-p_[i]) + cobs*path_l[i]*log(tau) - cobs*path_l[i]*log(1-tau+exp(-p_s))
        
      f_[i+1] <- lambda*score_[i] + delta*f_[i]
    }
    loglike <- sum(loglike)
    return(-loglike)
  }
  
  ###Estimation###
  
  #Starting Values
  
  param_0 <- c(f1 = 0 , w = 1, alpha = 0.04, beta = 0.3, lambda = 0.0002, delta = 0.05)  
  
  
  #COMPLETE CYCLE
  pe_com <- matrix(data = NA, nrow = simulations, ncol = length(param_0))
  tobs <- nrow(y_m)
  y_m <- as.matrix(y_m)
  
  for(i in 1:simulations){
    fit <- optim(par = param_0, fn = loglikelihood, method = "BFGS", path_l =  gen_success_rate_com[,i],
                 path_n = as.numeric(y_m[,i]), tobs = tobs, cobs = cobs)
    pe_com[i,] <- fit$par 
    print(i)
  }
  
  #Compute bias for each parameter
  #pe_com_bias[,1] <- pe_com[,1] - true_param[1]
  #pe_com_bias[,2] <- pe_com[,2] - true_param[2]
  #pe_com_bias[,3] <- pe_com[,3] - true_param[3]
  #pe_com_bias[,4] <- pe_com[,4] - true_param[4]
  #pe_com_bias[,5] <- pe_com[,5] - true_param[5]
  #pe_com_bias[,6] <- pe_com[,6] - true_param[6]
  
  
  #2C, 2R
  param_0 <- c(f1 = 0 , w = 1, alpha = 0.04, beta = 0.03, lambda = 0.00002, delta = 0.005) 
  pe_2c2r <- matrix(data = NA, nrow = simulations, ncol = length(param_0))
  tobs <- nrow(y_m_2c2r)
  y_m_2c2r <- as.matrix(y_m_2c2r)
  
  for(i in 1:simulations){
    fit <- optim(par = param_0, fn = loglikelihood, method = "BFGS", path_l =  gen_success_rate_2c2r[,i],
                 path_n = as.numeric(y_m_2c2r[,i]), tobs = tobs, cobs = cobs)
    pe_2c2r[i,] <- fit$par 
    print(i)
  }
  
  
  #2C, 2SG
  pe_2c2sg <- matrix(data = NA, nrow = simulations, ncol = length(param_0))
  tobs <- nrow(y_m_2c2sg)
  y_m_2c2sg <- as.matrix(y_m_2c2sg)
  
  for(i in 1:simulations){
    fit <- optim(par = param_0, fn = loglikelihood, method = "BFGS", path_l =  gen_success_rate_2c2sg[,i],
                 path_n = as.numeric(y_m_2c2sg[,i]), tobs = tobs, cobs = cobs)
    pe_2c2sg[i,] <- fit$par 
    print(i)
  }
  
  
  #2R, 2SG
  pe_2r2sg <- matrix(data = NA, nrow = simulations, ncol = length(param_0))
  tobs <- nrow(y_m_2r2sg)
  y_m_2r2sg <- as.matrix(y_m_2r2sg)
  
  for(i in 1:simulations){
    fit <- optim(par = param_0, fn = loglikelihood, method = "BFGS", path_l =  gen_success_rate_2r2sg[,i],
                 path_n = as.numeric(y_m_2r2sg[,i]), tobs = tobs, cobs = cobs)
    pe_2r2sg[i,] <- fit$par 
    print(i)
  }
  
  
  
  #2C
  pe_2c <- matrix(data = NA, nrow = simulations, ncol = length(param_0))
  tobs <- nrow(y_m_2c)
  y_m_2c <- as.matrix(y_m_2c)
  
  for(i in 1:simulations){
    fit <- optim(par = param_0, fn = loglikelihood, method = "BFGS", path_l =  gen_success_rate_2c[,i],
                 path_n = as.numeric(y_m_2c[,i]), tobs = tobs, cobs = cobs)
    pe_2c[i,] <- fit$par 
    print(i)
  }
  
  
  
  #2R
  pe_2r <- matrix(data = NA, nrow = simulations, ncol = length(param_0))
  tobs <- nrow(y_m_2r)
  y_m_2r <- as.matrix(y_m_2r)
  
  for(i in 1:simulations){
    fit <- optim(par = param_0, fn = loglikelihood, method = "BFGS", path_l =  gen_success_rate_2r[,i],
                 path_n = as.numeric(y_m_2r[,i]), tobs = tobs, cobs = cobs)
    pe_2r[i,] <- fit$par 
    print(i)
  }
  
  
  #2SG
  pe_2sg <- matrix(data = NA, nrow = simulations, ncol = length(param_0))
  tobs <- nrow(y_m_2sg)
  y_m_2sg <- as.matrix(y_m_2sg)
  
  for(i in 1:simulations){
    fit <- optim(par = param_0, fn = loglikelihood, method = "BFGS", path_l =  gen_success_rate_2sg[,i],
                 path_n = as.numeric(y_m_2sg[,i]), tobs = tobs, cobs = cobs)
    pe_2sg[i,] <- fit$par 
    print(i)
  }
  
  bias_com_GAS <- pe_com - true_param 
  bias_2r2sg_GAS <- pe_2r2sg - true_param 
  bias_2c2sg_GAS <- pe_2c2sg - true_param 
  bias_2c2r_GAS <- pe_2c2r - true_param 
  bias_2c_GAS <- pe_2c - true_param 
  bias_2r_GAS <- pe_2r - true_param 
  bias_2sg_GAS <- pe_2sg - true_param  
  
  
  #Tests for significant difference in bias -> Compare all to Complte sample
  t.test(bias_com_GAS[,1], bias_2r2sg_GAS[,1], alternative = "two.sided", conf.level = 0.95)
  t.test(bias_com_GAS[,1], bias_2c2sg_GAS[,1], alternative = "two.sided", conf.level = 0.95)
  t.test(bias_com_GAS[,1], bias_2c2r_GAS[,1], alternative = "two.sided", conf.level = 0.95)
  t.test(bias_com_GAS[,1], bias_2c_GAS[,1], alternative = "two.sided", conf.level = 0.95)
  t.test(bias_com_GAS[,1], bias_2r_GAS[,1], alternative = "two.sided", conf.level = 0.95)
  t.test(bias_com_GAS[,1], bias_2sg_GAS[,1], alternative = "two.sided", conf.level = 0.95)
  
  t.test(bias_com_GAS[,2], bias_2r2sg_GAS[,2], alternative = "two.sided", conf.level = 0.95)
  t.test(bias_com_GAS[,2], bias_2c2sg_GAS[,2], alternative = "two.sided", conf.level = 0.95)
  t.test(bias_com_GAS[,2], bias_2c2r_GAS[,2], alternative = "two.sided", conf.level = 0.95)
  t.test(bias_com_GAS[,2], bias_2c_GAS[,2], alternative = "two.sided", conf.level = 0.95)
  t.test(bias_com_GAS[,2], bias_2r_GAS[,2], alternative = "two.sided", conf.level = 0.95)
  t.test(bias_com_GAS[,2], bias_2sg_GAS[,2], alternative = "two.sided", conf.level = 0.95)
  
  t.test(bias_com_GAS[,3], bias_2r2sg_GAS[,3], alternative = "two.sided", conf.level = 0.95)
  t.test(bias_com_GAS[,3], bias_2c2sg_GAS[,3], alternative = "two.sided", conf.level = 0.95)
  t.test(bias_com_GAS[,3], bias_2c2r_GAS[,3], alternative = "two.sided", conf.level = 0.95)
  t.test(bias_com_GAS[,3], bias_2c_GAS[,3], alternative = "two.sided", conf.level = 0.95)
  t.test(bias_com_GAS[,3], bias_2r_GAS[,3], alternative = "two.sided", conf.level = 0.95)
  t.test(bias_com_GAS[,3], bias_2sg_GAS[,3], alternative = "two.sided", conf.level = 0.95)
  
  t.test(bias_com_GAS[,4], bias_2r2sg_GAS[,4], alternative = "two.sided", conf.level = 0.95)
  t.test(bias_com_GAS[,4], bias_2c2sg_GAS[,4], alternative = "two.sided", conf.level = 0.95)
  t.test(bias_com_GAS[,4], bias_2c2r_GAS[,4], alternative = "two.sided", conf.level = 0.95)
  t.test(bias_com_GAS[,4], bias_2c_GAS[,4], alternative = "two.sided", conf.level = 0.95)
  t.test(bias_com_GAS[,4], bias_2r_GAS[,4], alternative = "two.sided", conf.level = 0.95)
  t.test(bias_com_GAS[,4], bias_2sg_GAS[,4], alternative = "two.sided", conf.level = 0.95)
  
  t.test(bias_com_GAS[,5], bias_2r2sg_GAS[,5], alternative = "two.sided", conf.level = 0.95)
  t.test(bias_com_GAS[,5], bias_2c2sg_GAS[,5], alternative = "two.sided", conf.level = 0.95)
  t.test(bias_com_GAS[,5], bias_2c2r_GAS[,5], alternative = "two.sided", conf.level = 0.95)
  t.test(bias_com_GAS[,5], bias_2c_GAS[,5], alternative = "two.sided", conf.level = 0.95)
  t.test(bias_com_GAS[,5], bias_2r_GAS[,5], alternative = "two.sided", conf.level = 0.95)
  t.test(bias_com_GAS[,5], bias_2sg_GAS[,5], alternative = "two.sided", conf.level = 0.95)
  
  t.test(bias_com_GAS[,6], bias_2r2sg_GAS[,6], alternative = "two.sided", conf.level = 0.95)
  t.test(bias_com_GAS[,6], bias_2c2sg_GAS[,6], alternative = "two.sided", conf.level = 0.95)
  t.test(bias_com_GAS[,6], bias_2c2r_GAS[,6], alternative = "two.sided", conf.level = 0.95)
  t.test(bias_com_GAS[,6], bias_2c_GAS[,6], alternative = "two.sided", conf.level = 0.95)
  t.test(bias_com_GAS[,6], bias_2r_GAS[,6], alternative = "two.sided", conf.level = 0.95)
  t.test(bias_com_GAS[,6], bias_2sg_GAS[,6], alternative = "two.sided", conf.level = 0.95)
  
  
  
  #Quantigying  
  bias_2sg_GAS <- as.data.table(bias_2sg_GAS)
  bias_2sg_GAS[, SG := 1]
  bias_2sg_GAS[, C := 0]
  bias_2sg_GAS[, R := 0]
  bias_2sg_GAS[, TwoS := 0]
  bias_2sg_GAS[, Complete := 0]
  
  bias_2r_GAS <- as.data.table(bias_2r_GAS)
  bias_2r_GAS[, SG := 0]
  bias_2r_GAS[, C := 0]
  bias_2r_GAS[, R := 1]
  bias_2r_GAS[, TwoS := 0]
  bias_2r_GAS[, Complete := 0]
  
  bias_2c_GAS <- as.data.table(bias_2c_GAS)
  bias_2c_GAS[, SG := 0]
  bias_2c_GAS[, C := 1]
  bias_2c_GAS[, R := 0]
  bias_2c_GAS[, TwoS := 0]
  bias_2c_GAS[, Complete := 0]
  
  bias_2r2sg_GAS <- as.data.table(bias_2r2sg_GAS)
  bias_2r2sg_GAS[, SG := 1]
  bias_2r2sg_GAS[, C := 0]
  bias_2r2sg_GAS[, R := 1]
  bias_2r2sg_GAS[, TwoS := 1]
  bias_2r2sg_GAS[, Complete := 0]
  
  bias_2c2sg_GAS <- as.data.table(bias_2c2sg_GAS)
  bias_2c2sg_GAS[, SG := 1]
  bias_2c2sg_GAS[, C := 1]
  bias_2c2sg_GAS[, R := 0]
  bias_2c2sg_GAS[, TwoS := 1]
  bias_2c2sg_GAS[, Complete := 0]
  
  bias_2c2r_GAS <- as.data.table(bias_2c2r_GAS)
  bias_2c2r_GAS[, SG := 0]
  bias_2c2r_GAS[, C := 1]
  bias_2c2r_GAS[, R := 1]
  bias_2c2r_GAS[, TwoS := 1]
  bias_2c2r_GAS[, Complete := 0]
  
  bias_com_GAS <- as.data.table(bias_com_GAS)
  bias_com_GAS[, SG := 1]
  bias_com_GAS[, C := 1]
  bias_com_GAS[, R := 1]
  bias_com_GAS[, TwoS := 0]
  bias_com_GAS[, Complete := 0]
  
  bias_all_GAS <- rbind(bias_com_GAS, bias_2r2sg_GAS, bias_2c2r_GAS, bias_2c2sg_GAS, bias_2sg_GAS, bias_2r_GAS, bias_2c_GAS)
  quant_bias_GAS_1 <- lm(V1 ~ SG + C + R, data = bias_all_GAS)
  quant_bias_GAS_2 <- lm(V2 ~ SG + C + R, data = bias_all_GAS)
  quant_bias_GAS_3 <- lm(V3 ~ SG + C + R , data = bias_all_GAS)
  quant_bias_GAS_4 <- lm(V4 ~ SG + C + R, data = bias_all_GAS)
  quant_bias_GAS_5 <- lm(V5 ~ SG + C + R , data = bias_all_GAS)
  quant_bias_GAS_6 <- lm(V6 ~ SG + C + R , data = bias_all_GAS)
  
  #Previous to using this function save the data collected so far
  summary.lm(quant_bias_GAS_1)
  summary.lm(quant_bias_GAS_2)
  summary.lm(quant_bias_GAS_3)
  summary.lm(quant_bias_GAS_4)
  summary.lm(quant_bias_GAS_5)
  summary.lm(quant_bias_GAS_6)
  
  
  
  par(mfrow = c(1,2))
  simulations_com$V1 <- as.numeric(simulations_com$V1)
  simulations_com$V501 <- as.numeric(simulations_com$V501)
  gdp_p <- ggplot(simulations_com, aes(x = 1:120, y = V1, group = 1, colour = States)) + 
    theme(plot.title = element_text(hjust = 0.5)) + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))  +
    geom_line() + 
    xlab("Observation") + 
    ylab("Risk Factor Values")
  
  pd_p <- ggplot(simulations_com, aes(x = 1:120, y = V501, group = 1, colour = States)) + 
    theme(plot.title = element_text(hjust = 0.5)) + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))  +
    geom_line() + 
    xlab("Observation") + 
    ylab("Default Rate Values")
  
  multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
    library(grid)
    
    # Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)
    
    numPlots = length(plots)
    
    # If layout is NULL, then use 'cols' to determine layout
    if (is.null(layout)) {
      # Make the panel
      # ncol: Number of columns of plots
      # nrow: Number of rows needed, calculated from # of cols
      layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                       ncol = cols, nrow = ceiling(numPlots/cols))
    }
    
    if (numPlots==1) {
      print(plots[[1]])
      
    } else {
      # Set up the page
      grid.newpage()
      pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
      
      # Make each plot, in the correct location
      for (i in 1:numPlots) {
        # Get the i,j matrix positions of the regions that contain this subplot
        matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
        
        print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                        layout.pos.col = matchidx$col))
      }
    }
  }
  
  multiplot(gdp_p, pd_p)
  