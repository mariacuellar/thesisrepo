# Analyzing JPAL data
# 2017-10-05


# Packages
# install.packages("randomForest")
#install.packages("doParallel")
library(randomForest)
# library(doParallel)
# detectCores()
# registerDoParallel(4)


# Set working directory
mywd = "/Users/mariacuellar/Desktop/CMU/Papers/2nd Heinz paper:ADA/Shaken Baby Syndrome/THESIS/Drafts, proposal, defense, notes/Notes/Simulations/"
setwd(mywd)

# Loading data
dat.stove = read.table("Stovedata.csv", head=TRUE, sep=",")

# Making Y, A, X1, X2 for code

names(dat.stove)
dat.stove$Y = dat.stove$cold_or_cough_BL
dat.stove$Y = as.numeric(dat.stove$Y)
dat.stove$Y = ifelse(dat.stove$Y==2, 0, dat.stove$Y)
dat.stove$Y = ifelse(dat.stove$Y==1, NA, dat.stove$Y)
dat.stove$Y = ifelse(dat.stove$Y==3, 1, dat.stove$Y)
dat.stove$A = dat.stove$treatment
dat.stove$X1 = dat.stove$male
dat.stove$X2 = dat.stove$primarycook
head(dat.stove)





#  Influence-function-based estimator, estimate nuisance parameters with plugin nonparametric untransformed Xs ----
dff = dat.stove

#  Defining the nuisance parameters
IF_N_X_mu0 = randomForest(Y~X1+X2, data = dff[which(dff$A==0),], na.action=na.exclude, type=regression)
IF_N_X_mu1 = randomForest(Y~X1+X2, data = dff[which(dff$A==1),], na.action=na.exclude, type=regression)
IF_N_X_pi = randomForest(A~X1+X2, data = dff, na.action=na.exclude, type=regression)

# Getting fitted values
IF_N_X_mu0_hat = as.numeric(predict(IF_N_X_mu0, newdata=dff))
IF_N_X_mu1_hat = as.numeric(predict(IF_N_X_mu1, newdata=dff))
IF_N_X_pi_hat = as.numeric(predict(IF_N_X_pi, newdata=dff))

# #  Defining the nuisance parameters
# IF_N_X_mu0 = npreg(Y~X1+X2, data = dff[which(dff$A==0),], tol=0.1, ftol=0.1)
# IF_N_X_mu1 = npreg(Y~X1+X2, data = dff[which(dff$A==1),], tol=0.1, ftol=0.1)
# IF_N_X_pi = npreg(A~X1+X2, data = dff, tol=0.1, ftol=0.1)
# 
# # Getting fitted values
# IF_N_X_mu0_hat = predict(IF_N_X_mu0, newdata = dff)
# IF_N_X_mu1_hat = predict(IF_N_X_mu1, newdata = dff)
# IF_N_X_pi_hat = predict(IF_N_X_pi, newdata = dff)
# plot(IF_N_X_mu0_hat, dff$Y)

# Defining my pseudo-outcome
IF_N_X_y = (1/IF_N_X_mu1_hat)*((IF_N_X_mu0_hat/IF_N_X_mu1_hat)*(1/IF_N_X_pi_hat)*A*(Y-IF_N_X_mu1_hat) - 
                                 (1/(1-IF_N_X_pi_hat))*(1-A)*(Y-IF_N_X_mu0_hat)) + (IF_N_X_mu1_hat-IF_N_X_mu0_hat)/IF_N_X_mu1_hat

# Fitting model
IF_N_X_model = try(nls(
  IF_N_X_y ~ fun.expit(beta1*X1 + beta2*X2),
  start=list(beta1=0, beta2=0),
  lower=c(-2, -2, -2, -2),
  upper=c(2, 2, 2, 2),
  algorithm="port",
  data=dff, nls.control(maxiter = 500))
  , silent=TRUE)

# Getting predicted values
IF_N_X_gammahat = try(predict(IF_N_X_model), silent=TRUE)
IF_N_X_gammahat_mean = mean(sample(IF_N_X_gammahat, size = 100))

# Coverage (Estimate coverage using the bootstrap)
# Function to estimate gamma
fun.IF_N_X_boot = function(xthedata){
  #  Defining the nuisance parameters
  IF_N_X_mu0 = randomForest(Y~X1+X2, data = xthedata[which(xthedata$A==0),], na.action=na.exclude, type=regression)
  IF_N_X_mu1 = randomForest(Y~X1+X2, data = xthedata[which(xthedata$A==1),], na.action=na.exclude, type=regression)
  IF_N_X_pi = randomForest(A~X1+X2, data = xthedata, na.action=na.exclude, type=regression)
  
  # Getting fitted values
  IF_N_X_mu0_hat = as.numeric(predict(IF_N_X_mu0, newdata=xthedata))
  IF_N_X_mu1_hat = as.numeric(predict(IF_N_X_mu1, newdata=xthedata))
  IF_N_X_pi_hat = as.numeric(predict(IF_N_X_pi, newdata=xthedata))
  
  # Defining my pseudo-outcome
  IF_N_X_y = (1/IF_N_X_mu1_hat)*((IF_N_X_mu0_hat/IF_N_X_mu1_hat)*(1/IF_N_X_pi_hat)*A*(Y-IF_N_X_mu1_hat) - 
                                   (1/(1-IF_N_X_pi_hat))*(1-A)*(Y-IF_N_X_mu0_hat)) + (IF_N_X_mu1_hat-IF_N_X_mu0_hat)/IF_N_X_mu1_hat
  
  # Fitting model
  IF_N_X_model = try(nls(
    IF_N_X_y ~ fun.expit(beta1*X1 + beta2*X2),
    start=list(beta1=0, beta2=0),
    lower=c(-2, -2, -2, -2),
    upper=c(2, 2, 2, 2),
    algorithm="port",
    data=xthedata, nls.control(maxiter = 500))
    , silent=TRUE)
  
  # Getting predicted values
  IF_N_X_gammahat = try(fun.expit(predict(IF_N_X_model)), silent=TRUE)
  IF_N_X_gammahat_mean = mean(sample(IF_N_X_gammahat, size = 100))
  
  return(IF_N_X_gammahat_mean)
}

# Bootstrap
bootreps1=10
IF_N_X_bootvec = 0
#stime = system.time({ # this is for testing the time. It's 3-4 times faster with dopar than without.
IF_N_X_bootvec.1 = foreach(i=1:bootreps1) %dopar% {
  index.b = sample(x = 1:nrow(dff), replace=TRUE) # randomize indices
  newdff = dff[index.b,] # select new sample of rows from dff
  IF_N_X_gammahat_star = fun.IF_N_X_boot(newdff) # estimating gammahat for new dataset
  IF_N_X_bootvec[i] <- IF_N_X_gammahat_star
} 
#})
#stime

# Standard error
IF_N_X_gammahat_sd = sd(unlist(IF_N_X_bootvec.1))

# Confidence interval
IF_N_X_ci_lb = IF_N_X_gammahat_mean - 2*IF_N_X_gammahat_sd / sqrt(samplesize)
IF_N_X_ci_ub = IF_N_X_gammahat_mean + 2*IF_N_X_gammahat_sd / sqrt(samplesize)
IF_N_X_ci = paste(IF_N_X_ci_lb, IF_N_X_ci_ub, sep=", ")

# "0.223827896627772, 0.224211845533956"







# Plug-in, parametric, transformed X's (misspecified model)

# defining the nuisance parameters
PI_P_X_mu0 = glm(Y~X1+X2, data = dff[which(dff$A==0),], family = "binomial")
PI_P_X_mu1 = glm(Y~X1+X2, data = dff[which(dff$A==1),], family = "binomial")

# getting fitted values
PI_P_X_mu0_hat = predict(PI_P_X_mu0, newdata = dff, type="response")
PI_P_X_mu1_hat = predict(PI_P_X_mu1, newdata = dff, type="response")

# calculating gamma hat
PI_P_X_gammahat = (PI_P_X_mu1_hat - PI_P_X_mu0_hat)/PI_P_X_mu1_hat
PI_P_X_gammahat_mean = mean(as.numeric(sample(PI_P_X_gammahat, size = 100)), na.rm=TRUE)

# Coverage (Estimate coverage using the bootstrap)
# Function to estimate gamma
xthedata=dff
fun.PI_P_X_boot = function(xthedata){
  # defining the nuisance parameters
  PI_P_X_mu0 = glm(Y~X1+X2, data = xthedata[which(xthedata$A==0),], family = "binomial")
  PI_P_X_mu1 = glm(Y~X1+X2, data = xthedata[which(xthedata$A==1),], family = "binomial")
  
  # getting fitted values
  PI_P_X_mu0_hat = predict(PI_P_X_mu0, newdata = xthedata, type="response")
  PI_P_X_mu1_hat = predict(PI_P_X_mu1, newdata = xthedata, type="response")
  
  # calculating gamma hat
  PI_P_X_gammahat = (PI_P_X_mu1_hat - PI_P_X_mu0_hat)/PI_P_X_mu1_hat
  PI_P_X_gammahat_mean = mean(as.numeric(sample(PI_P_X_gammahat, size = 100)), na.rm=TRUE)
  
  return(PI_P_X_gammahat_mean)
}

# Bootstrap 
bootreps=100
PI_P_X_bootvec <- 0
PI_P_X_bootvec.1 = foreach(i=1:bootreps) %dopar% {
  index.b = sample(x = 1:nrow(dff), replace=TRUE) # randomize indices
  newdff = dff[index.b,] # select new sample of rows from dff
  PI_P_X_gammahat_star = fun.PI_P_X_boot(newdff) # estimating gammahat for new dataset
  PI_P_X_bootvec[i] <- PI_P_X_gammahat_star
}

# Standard error
PI_P_X_gammahat_sd = sd(unlist(PI_P_X_bootvec.1))

# Confidence interval
PI_P_X_ci_lb = PI_P_X_gammahat_mean - 2*PI_P_X_gammahat_sd / sqrt(samplesize)
PI_P_X_ci_ub = PI_P_X_gammahat_mean + 2*PI_P_X_gammahat_sd / sqrt(samplesize)
PI_P_X_ci = paste(PI_P_X_ci_lb, PI_P_X_ci_ub, sep=", ")

"0.0183518423488624, 0.0199616313618164"






