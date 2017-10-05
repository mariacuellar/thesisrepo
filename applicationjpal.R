# # Thesis simulations, doc for Cosma
# # October 2, 2017
# # Making my own simulation for the parameter of interest: PC = gamma = 1-1/RR = 1-mu0/mu1.


# Packages

# install.packages("np")
# install.packages("glmnet")
# install.packages("boot")
# install.packages("parallel")
# install.packages("doParallel")

library(np)
library(glmnet)
library(boot)
library(parallel)
library(doParallel)



# Definte expit and logit functions
fun.expit = function(x){ return( exp(x)/(1+exp(x)) ) }
fun.logit = function(x){ return( log(x)/log(1-x) ) }


# Function to simulate estimation of PC (gamma) as plugin (parametric and nonparametric), 
# and as an IF where the nuisance parameters are parametric
samplesize=2000
fun.datageneration.with4 = function(samplesize){
  # let us know how far it's gone
  cat(paste(samplesize,"... "))
  
  # true parameters
  index = 1:samplesize
  beta = 0.5
  Y1 = rbinom(n = samplesize, size = 1, prob = beta)
  
  # correct model
  X1 = rnorm(n = samplesize, mean = 0, sd = 1) 
  X2 = rnorm(n = samplesize, mean = 0, sd = 1)
  X3 = rnorm(n = samplesize, mean = 0, sd = 1)
  X4 = rnorm(n = samplesize, mean = 0, sd = 1)
  
  # misspecified model
  X1star = exp(X1/2) 
  X2star = X2/(1 + exp(X1)) + 10
  X3star = (X1*X3/25 + 0.6)^3
  X4star = (X2 + X4 + 20)^2
  
  pi = fun.expit(-X1+0.5*X2-0.25*X3-0.1*X4)
  
  A = rbinom(n = samplesize, size = 1, prob = fun.expit(-X1+0.5*X2-0.25*X3-0.1*X4))
  
  # Defining my parameter
  beta = 0.5
  mu0 = beta/(1+exp(-X1+0.5*X2-0.25*X3-0.1*X4))
  mu1 = rep(beta, samplesize)
  gamma = fun.expit(-X1+0.5*X2-0.25*X3-0.1*X4) #gamma = 1 - mu0/mu1
  meangamma = mean(sample(gamma, 100))
  
  # generate data frame
  df = as.data.frame(cbind(index, X1, X2, X3, X4, X1star, X2star, X3star, X4star, pi, gamma, A, Y1))
  head(df)
  
  # generate Y0 conditional on combinations of values of A and Y1
  dfy11 = df[which(df$Y1==1),]
  dfy11$Y0 = rbinom(n = nrow(dfy11) , size = 1, prob = (1-gamma)) # or is it location = fun.expit(t(PSI)*X), scale = 0?
  
  dfy10 = df[which(df$Y1==0),]
  dfy10$Y0 = 0
  
  # add Y0 to dataframe
  df_wy0 = as.data.frame(rbind(dfy11, dfy10))
  
  # apply consistency to get Y
  df_wy0$Y = ifelse(df_wy0$A==1, df_wy0$Y1, df_wy0$Y0) 
  Y = df_wy0$Y
  
  # ordering data so it's as it was at the beginning
  dff = df_wy0[order(df_wy0$index),] 
  
  return(dff)
}


# Fitting it with one parameter
fun.obj = function(param, theY){
  #betapointfive = exp(param[1])
  beta1 = param[1]
  Sel = A==0
  phat = 0.5/(1+exp(beta1*X1[Sel]))
  dens = dbinom(Y[Sel], 1, phat, log=TRUE)
  return(-sum(dens))
}

fun.obj(c(log(0.5),0), Y)


result = nlm(fun.obj, 0, Y, print.level=2, iterlim=10)
predict(result)

estimate = result[2][[1]]

.5/(1+exp(estimate))














fun.datageneration.with2 = function(samplesize){
  # let us know how far it's gone
  cat(paste(samplesize,"... "))
  
  # true parameters
  index = 1:samplesize
  beta = 0.5
  Y1 = rbinom(n = samplesize, size = 1, prob = beta)
  
  # correct model
  X1 = rnorm(n = samplesize, mean = 0, sd = 1) 
  X2 = rnorm(n = samplesize, mean = 0, sd = 1)
  
  # misspecified model
  X1star = exp(X1/2) 
  X2star = X2/(1 + exp(X1)) + 10
  
  pi = fun.expit(-X1+0.5*X2)
  
  A = rbinom(n = samplesize, size = 1, prob = fun.expit(-X1+0.5*X2))
  
  # Defining my parameter
  beta = 0.5
  mu0 = 0.5/(1+exp(-X1+0.5*X2))
  mu1 = rep(beta, samplesize)
  gamma = fun.expit(-X1+0.5*X2) #gamma = 1 - mu0/mu1
  meangamma = mean(sample(gamma, 100))
  
  # generate data frame
  df = as.data.frame(cbind(index, X1, X2, X1star, X2star, pi, gamma, A, Y1))
  head(df)
  
  # generate Y0 conditional on combinations of values of A and Y1
  dfy11 = df[which(df$Y1==1),]
  dfy11$Y0 = rbinom(n = nrow(dfy11) , size = 1, prob = (1-gamma)) # or is it location = fun.expit(t(PSI)*X), scale = 0?
  
  dfy10 = df[which(df$Y1==0),]
  dfy10$Y0 = 0
  
  # add Y0 to dataframe
  df_wy0 = as.data.frame(rbind(dfy11, dfy10))
  
  # apply consistency to get Y
  df_wy0$Y = ifelse(df_wy0$A==1, df_wy0$Y1, df_wy0$Y0) 
  Y = df_wy0$Y
  
  # ordering data so it's as it was at the beginning
  dff = df_wy0[order(df_wy0$index),] 
  
  return(dff)
}


dat2 = fun.datageneration.with2(2000)
attach(dat2)

# fitting it with two parameters
fun.obj = function(param, theY){
  #betapointfive = exp(param[1])
  beta1 = param[1]
  beta2 = param[2]
  Sel = A==0
  phat = 0.5/(1+exp(beta1*X1[Sel] + beta2*X2[Sel]))
  dens = dbinom(Y[Sel], 1, phat, log=TRUE)
  return(-sum(dens))
}

fun.obj(c(log(0.5),0), Y)


result = nlm(fun.obj, c(0,0), Y, print.level=2, iterlim=10)
result[2]

.5/(1+exp(estimate))


install.packages("devtools")
library(devtools)
x=c(1,2,3,4)
dx.expit(x)
??dx.expit

dat = fun.datageneration.with2(2000)

# Wrong way to fit binary
PI_P_X_mu0 = nls(Y~.5/(1+exp(beta1*X1+beta2*X2)), 
                 start=list(beta1=-1, beta2=-.25), 
                 data = dat[which(dat$A==0),])

PI_P_X_mu0_hat = predict(PI_P_X_mu0, newdata = dat)
plot(mu0, PI_P_X_mu0_hat)

PI_P_X_mu1 = glm(Y~X1+X2+X3+X4, data = dff[which(dff$A==1),], family = "binomial")

# getting fitted values (after inverse link function)
PI_P_X_mu0_hat = predict(PI_P_X_mu0, newdata = dff)
PI_P_X_mu1_hat = predict(PI_P_X_mu1, newdata = dff, type="response")

# plot(PI_P_X_mu0_hat, mu0) # Why are my predictions so bad??

# calculating gamma hat
PI_P_X_gammahat = (PI_P_X_mu1_hat - PI_P_X_mu0_hat)/PI_P_X_mu1_hat
PI_P_X_gammahat_mean = mean(sample(PI_P_X_gammahat, size = 100))

# RMSE
PI_P_X_RMSE = sqrt(  mean( (PI_P_X_gammahat_mean - meangamma)^2 )  ) #

# Bias
PI_P_X_bias = mean(abs(PI_P_X_gammahat_mean - meangamma))

# Coverage (Estimate coverage using the bootstrap)
# Function to estimate gamma
fun.PI_P_X_boot = function(xthedata){
  # Plug-in, parametric, untransformed X's (correctly specified model)
  PI_P_X_mu0 = nls(Y~beta1/(1+exp(-X1+0.5*X2-0.25*X3-0.1*X4)), start=list(beta1=0.5), data = xthedata[which(xthedata$A==0),]) # take away the parameters and let it estimate them
  PI_P_X_mu1 = glm(Y~X1+X2+X3+X4, data = xthedata[which(xthedata$A==1),], family = "binomial")
  
  # getting fitted values (after inverse link function)
  PI_P_X_mu0_hat = predict(PI_P_X_mu0, newdata = xthedata)
  PI_P_X_mu1_hat = predict(PI_P_X_mu1, newdata = xthedata, type="response")
  
  # calculating gamma hat
  PI_P_X_gammahat = (PI_P_X_mu1_hat - PI_P_X_mu0_hat)/PI_P_X_mu1_hat
  return(PI_P_X_gammahat)
}


# Bootstrap
bootreps=1000
PI_P_X_bootvec = 0
#stime = system.time({ # this is for testing the time. It's 3 times faster with dopar than without.
PI_P_X_bootvec.1 = foreach(i=1:bootreps, .options.multicore=list(preschedule=TRUE)) %dopar% {
  index.b = sample(x = 1:nrow(dff), replace=FALSE) # randomize indices
  newdff = dff[index.b,] # select new sample of rows from dff
  PI_P_X_gammahat_star = fun.PI_P_X_boot(newdff) # estimating gammahat for new dataset
  PI_P_X_gammahat_star_mean = mean(sample(PI_P_X_gammahat_star, size = 100)) # getting a sample of 100 gamma hat stars and taking their average
  PI_P_X_bootvec[i] <- PI_P_X_gammahat_star_mean
}
#})
#stime

# Standard error
PI_P_X_gammahat_sd = sd(unlist(PI_P_X_bootvec.1))

# Confidence interval
PI_P_X_ci_lb = PI_P_X_gammahat_mean - 2*PI_P_X_gammahat_sd / sqrt(samplesize)
PI_P_X_ci_ub = PI_P_X_gammahat_mean + 2*PI_P_X_gammahat_sd / sqrt(samplesize)
PI_P_X_ci = paste(PI_P_X_ci_lb, PI_P_X_ci_ub, sep=", ")












# 2) Misspecified plug-in parametric model
PI_P_Xstar_mu0 = glm(Y~X1star+X2star+X3star+X4star, data = dff[which(dff$A==0),], family = "binomial")
PI_P_Xstar_mu1 = glm(Y~X1star+X2star+X3star+X4star, data = dff[which(dff$A==1),], family = "binomial")

# getting fitted values (after inverse link function)
PI_P_Xstar_mu0_hat = predict(PI_P_Xstar_mu0, newdata = dff, type="response")
PI_P_Xstar_mu1_hat = predict(PI_P_Xstar_mu1, newdata = dff, type="response")

# calculating gamma hat
PI_P_Xstar_gammahat  = (PI_P_Xstar_mu1_hat - PI_P_Xstar_mu0_hat)/PI_P_Xstar_mu1_hat
PI_P_Xstar_gammahat_mean = mean(sample(PI_P_Xstar_gammahat, size = 100))

# RMSE 
PI_P_Xstar_RMSE = sqrt(  mean( (PI_P_Xstar_gammahat_mean - meangamma)^2 )  )

# Bias
PI_P_Xstar_bias = mean(abs(PI_P_Xstar_gammahat_mean - meangamma))

# Coverage

# Coverage (Estimate coverage using the bootstrap)
# Function to estimate gamma
fun.PI_P_Xstar_boot = function(xthedata){
  PI_P_Xstar_mu0_hat = predict(glm(Y~X1star+X2star+X3star+X4star, data = xthedata[which(xthedata$A==0),], family = "binomial"), newdata = xthedata, type="response")
  PI_P_Xstar_mu1_hat = predict(glm(Y~X1star+X2star+X3star+X4star, data = xthedata[which(xthedata$A==1),], family = "binomial"), newdata = xthedata, type="response")
  PI_P_Xstar_gammahat = (PI_P_Xstar_mu1_hat - PI_P_Xstar_mu0_hat)/PI_P_Xstar_mu1_hat
  PI_P_Xstar_gammahat_mean = mean(sample(PI_P_Xstar_gammahat, size = 100))
  return(PI_P_Xstar_gammahat_mean)
}

# Bootstrap 
# bootreps defined above
PI_P_Xstar_bootvec <- 0
PI_P_Xstar_bootvec.1 = foreach(i=1:bootreps) %dopar% {
  index.b = sample(x = 1:nrow(dff), replace=TRUE) # randomize indices
  newdff = dff[index.b,] # select new sample of rows from dff
  PI_P_Xstar_gammahat_star = fun.PI_P_Xstar_boot(newdff) # estimating gammahat for new dataset
  PI_P_Xstar_bootvec[i] <- PI_P_Xstar_gammahat_star
}

# Standard error
PI_P_Xstar_gammahat_sd = sd(unlist(PI_P_Xstar_bootvec.1))

# Confidence interval
PI_P_Xstar_ci_lb = PI_P_Xstar_gammahat_mean - 2*PI_P_Xstar_gammahat_sd / sqrt(samplesize)
PI_P_Xstar_ci_ub = PI_P_Xstar_gammahat_mean + 2*PI_P_Xstar_gammahat_sd / sqrt(samplesize)
PI_P_Xstar_ci = paste(PI_P_Xstar_ci_lb, PI_P_Xstar_ci_ub, sep=", ")













# 3) Plug-in, nonparametric, untransformed X's estimator:

# Kernel----
# Fitting model, a bit faster than without tol, ftol, but still slow
# got bandwidths from running it once on a sample of the data
PI_N_X_mu0 = npreg(Y~X1+X2+X3+X4, data = dff[which(dff$A==0),], tol=0.1, ftol=0.1, 
                   bws=c(4760221, 403285.9, 1.333527, 2.563413))

PI_N_X_mu1 = npreg(Y~X1+X2+X3+X4, data = dff[which(dff$A==1),], tol=0.1, ftol=0.1,
                   bws=c(809346.9, 1743946, 0.0643948029, 806061.1))

# Getting fitted values
PI_N_X_mu0_hat = predict(PI_N_X_mu0, newdata = dff)
PI_N_X_mu1_hat = predict(PI_N_X_mu1, newdata = dff)

# plot(PI_N_X_mu0_hat, mu0) # Why are my predictions so bad??

# Gamma hat
PI_N_X_gammahat = (PI_N_X_mu1_hat - PI_N_X_mu0_hat)/PI_N_X_mu1_hat
PI_N_X_gammahat_mean = mean(sample(PI_N_X_gammahat, size = 100))

# RMSE 
PI_N_X_RMSE = sqrt(  mean( (PI_N_X_gammahat_mean - meangamma)^2 )  )

# Bias
PI_N_X_bias = mean(abs(PI_N_X_gammahat_mean - meangamma))

# Confidence interval
# # No valid confidence interval here.
PI_N_X_ci = NA

# # Coverage
# # No valid confidence interval here.
# PI_N_X_coverage = NA












# 4) Plug-in, nonparametric, transformed X's estimator:

# Kernel----
# Fitting model, a bit faster than without tol, ftol, but still slow
# got bandwidths from running it once on a sample of the data
PI_N_Xstar_mu0 = npreg(Y~X1star+X2star+X3star+X4star, data = dff[which(dff$A==0),], tol=0.1, ftol=0.1,
                       bws=c(5229101.686, 590245.5023, 0.02962718249, 46.19266321))

PI_N_Xstar_mu1 = npreg(Y~X1star+X2star+X3star+X4star, data = dff[which(dff$A==1),], tol=0.1, ftol=0.1,
                       bws=c(0.8807241131, 3242719.391, 0.0643948029, 85037300.36))

# Getting fitted values
PI_N_Xstar_mu0_hat = predict(PI_N_Xstar_mu0, newdata = dff)#$mean # with faster kernel
PI_N_Xstar_mu1_hat = predict(PI_N_Xstar_mu1, newdata = dff)#$mean

# Gamma hat
PI_N_Xstar_gammahat = (PI_N_Xstar_mu1_hat - PI_N_Xstar_mu0_hat)/PI_N_Xstar_mu1_hat
PI_N_Xstar_gammahat_mean = mean(sample(PI_N_Xstar_gammahat, size = 100))

# RMSE 
PI_N_Xstar_RMSE = sqrt(  mean( (PI_N_Xstar_gammahat_mean - meangamma)^2 )  )

# Bias
PI_N_Xstar_bias = mean(abs(PI_N_Xstar_gammahat_mean - meangamma))

# Confidence interval
# # No valid confidence interval here.
PI_N_Xstar_ci = NA

# # Coverage
# # No valid confidence interval here.
# PI_N_Xstar_coverage = NA















#  5) Influence-function-based estimator, estimate nuisance parameters with plugin parametric untransformed Xs ----
#  Defining the nuisance parameters
IF_P_X_mu0 = nls(Y~beta1/(1+exp(-X1+0.5*X2-0.25*X3-0.1*X4)), start=list(beta1=0.5), data = dff[which(dff$A==0),]) # same as PI_P_X_mu0_hat
IF_P_X_mu1 = glm(Y~X1+X2+X3+X4, data = dff[which(dff$A==1),], family = "binomial") # same as PI_P_X_mu1_hat
IF_P_X_pi = glm(A~X1+X2+X3+X4, data = dff, family = "binomial") # same as PI_P_X_mu1_hat

# getting fitted values (after inverse link function)
IF_P_X_mu0_hat = predict(IF_P_X_mu0, newdata = dff) # m0 # change this if you want to cheat by using true values!
IF_P_X_mu1_hat = predict(IF_P_X_mu1, newdata = dff, type="response") # mu1
IF_P_X_pi_hat = predict(IF_P_X_pi, newdata = dff, type="response") # pi

# Defining my pseudo-outcome
IF_P_X_ystar = (1/IF_P_X_mu1_hat)*((IF_P_X_mu0_hat/IF_P_X_mu1_hat)*(1/IF_P_X_pi_hat)*A*(Y-IF_P_X_mu1_hat) - 
                                     (1/(1-IF_P_X_pi_hat))*(1-A)*(Y-IF_P_X_mu0_hat)) + (IF_P_X_mu1_hat-IF_P_X_mu0_hat)/IF_P_X_mu1_hat

# Fitting model
IF_P_X_model = try(nls(
  IF_P_X_ystar ~ fun.expit(beta1*X1 + beta2*X2 + beta3*X3 + beta4*X4),
  start=list(beta1=0, beta2=0, beta3=0, beta4=0),
  lower=c(-2, -2, -2, -2), upper=c(2, 2, 2, 2),
  algorithm="port",
  data=dff, nls.control(maxiter = 500))
  , silent=TRUE)

# Getting predicted values
IF_P_X_gammahat = try(fun.expit(predict(IF_P_X_model)), silent=TRUE)
IF_P_X_gammahat_mean = mean(sample(IF_P_X_gammahat, size = 100))

# RMSE
IF_P_X_RMSE = try(sqrt(  mean( (IF_P_X_gammahat_mean - meangamma)^2 )  ), silent=TRUE)

# Show NA if I get an error message
IF_P_X_RMSE = ifelse(class(IF_P_X_RMSE)=="numeric", IF_P_X_RMSE, NA)

# Bias
IF_P_X_bias = try( mean(abs(IF_P_X_gammahat_mean - meangamma)) , silent=TRUE)

# Show NA if I get an error message
IF_P_X_bias = ifelse(class(IF_P_X_bias)=="numeric" && IF_P_X_bias!=Inf, IF_P_X_bias, NA)

# Coverage (Estimate coverage using the bootstrap)
# Function to estimate gamma
fun.IF_P_X_boot = function(xthedata){
  #  Defining the nuisance parameters
  IF_P_X_mu0 = nls(Y~beta1/(1+exp(-X1+0.5*X2-0.25*X3-0.1*X4)), start=list(beta1=0.5), data = xthedata[which(xthedata$A==0),]) # same as PI_P_X_mu0_hat
  IF_P_X_mu1 = glm(Y~X1+X2+X3+X4, data = xthedata[which(xthedata$A==1),], family = "binomial") # same as PI_P_X_mu1_hat
  IF_P_X_pi = glm(A~X1+X2+X3+X4, data = xthedata, family = "binomial") # same as PI_P_X_mu1_hat
  
  # getting fitted values (after inverse link function)
  IF_P_X_mu0_hat = predict(IF_P_X_mu0, newdata = dff) # m0 # change this if you want to cheat by using true values!
  IF_P_X_mu1_hat = predict(IF_P_X_mu1, newdata = dff, type="response") # mu1
  IF_P_X_pi_hat = predict(IF_P_X_pi, newdata = dff, type="response") # pi
  
  # Defining my pseudo-outcome
  IF_P_X_ystar = (1/IF_P_X_mu1_hat)*((IF_P_X_mu0_hat/IF_P_X_mu1_hat)*(1/IF_P_X_pi_hat)*A*(Y-IF_P_X_mu1_hat) - 
                                       (1/(1-IF_P_X_pi_hat))*(1-A)*(Y-IF_P_X_mu0_hat)) + (IF_P_X_mu1_hat-IF_P_X_mu0_hat)/IF_P_X_mu1_hat
  
  # Fitting model
  IF_P_X_model = try(nls(
    IF_P_X_ystar ~ fun.expit(beta1*X1 + beta2*X2 + beta3*X3 + beta4*X4),
    start=list(beta1=0, beta2=0, beta3=0, beta4=0),
    lower=c(-2, -2, -2, -2),
    upper=c(2, 2, 2, 2),
    algorithm="port",
    data=xthedata, nls.control(maxiter = 500))
    , silent=TRUE)
  
  # Getting predicted values for gamma hat
  IF_P_X_gammahat = try(fun.expit(predict(IF_P_X_model)), silent=TRUE)
  IF_P_X_gammahat_mean = mean(sample(IF_P_X_gammahat, size = 100))
  
  return(IF_P_X_gammahat_mean)
}

# Bootstrap
bootreps=100
IF_P_X_bootvec = 0
i=1
IF_P_X_bootvec.1 = foreach(i=1:bootreps) %dopar% {
  index.b = sample(x = 1:nrow(dff), replace=TRUE) # randomize indices
  newdff = dff[index.b,] # select new sample of rows from dff
  IF_P_X_gammahat_star = fun.IF_P_X_boot(newdff) # estimating gammahat for new dataset
  IF_P_X_bootvec[i] <- IF_P_X_gammahat_star
}

# Standard error
IF_P_X_gammahat_sd = sd(unlist(IF_P_X_bootvec.1))

# Confidence interval
IF_P_X_ci_lb = IF_P_X_gammahat_mean - 2*IF_P_X_gammahat_sd / sqrt(samplesize)
IF_P_X_ci_ub = IF_P_X_gammahat_mean + 2*IF_P_X_gammahat_sd / sqrt(samplesize)
IF_P_X_ci = paste(IF_P_X_ci_lb, IF_P_X_ci_ub, sep=", ")



# Function will return this
toreturn = c(PI_P_X_RMSE, PI_P_Xstar_RMSE, PI_N_X_RMSE, PI_N_Xstar_RMSE, IF_P_X_RMSE,
             PI_P_X_bias, PI_P_Xstar_bias, PI_N_X_bias, PI_N_Xstar_bias, IF_P_X_bias, 
             PI_P_X_ci, PI_P_Xstar_ci, PI_N_X_ci, PI_N_Xstar_ci, IF_P_X_ci
)

return(toreturn)
}



# this is not working right now, and I think the nls functions are to blame.
fun.simulate(2000) 

