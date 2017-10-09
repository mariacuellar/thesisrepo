# # Thesis simulations
# # Oct 8, 2017
# # Making my own simulation for the parameter of interest: PC = gamma = 1-1/RR = 1-mu0/mu1.

# # define working directory
WD_figs = "/home/mcuellar"
WD_thesis = "/home/mcuellar/Thesis"
WD_simulation = "/home/mcuellar/Thesis"

# Turn off warnings
#options(warn=-1)

# Turn on warnings
#options(warn=0)


# How many cores to use? # our server has 32, but James said I could only use 8.
# install.packages("doParallel")
# library(doParallel)
registerDoParallel(8)

# Definte expit and logit functions
fun.expit = function(x){ return( exp(x)/(1+exp(x)) ) }
fun.logit = function(x){ return( log(x)/log(1-x) ) }


# Function to simulate estimation of PC (gamma) as plugin, nonparametric, parametric IF, and nonparametric IF

samplesize = 1000 # 

# samplesize must be > 100.

# Data generation
fun.generate.data = function(samplesize){
  # set.seed(400) # seed for random number generator
  # let us know how far it's gone
  #cat(paste(samplesize,"... "))
  
  # true parameters
  index = 1:samplesize
  beta = 0.5
  Y1 = rbinom(n = samplesize, size = 1, prob = beta)
  
  X1 = rnorm(n = samplesize, mean = 0, sd = .001) # correct model
  X2 = rnorm(n = samplesize, mean = 0, sd = .001)
  X3 = rnorm(n = samplesize, mean = 0, sd = .001)
  X4 = rnorm(n = samplesize, mean = 0, sd = .001)
  
  X1star = exp(X1/2) # misspecified model
  X2star = X2/(1 + exp(X1)) + 10
  X3star = (X1*X3/25 + 0.6)^3
  X4star = (X2 + X4 + 20)^2
  
  pi = fun.expit(-X1+0.5*X2-0.25*X3-0.1*X4)
  
  A = rbinom(n = samplesize, size = 1, prob = fun.expit(-X1+0.5*X2-0.25*X3-0.1*X4))
  
  # Defining my parameter
  beta = 0.5
  mu0 = beta/(1+exp(-X1+0.5*X2-0.25*X3-0.1*X4))
  mu1 = rep(beta, samplesize)
  #gamma = 1 - mu0/mu1
  gamma = fun.expit(-X1+0.5*X2-0.25*X3-0.1*X4) # Don't we need to make it so that gamma is equal to 1-mu0/mu1?? It is.
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


  # 1. Plug-in, parametric, untransformed X's (correctly specified model)
fun.PI_P_X = function(dff){
  
  # data splitting
  index = sample(1:nrow(dff)) # randomize indices
  dff.train = dff[index,][1:(nrow(dff)/2),] # make training set
  dff.test = dff[index,][((nrow(dff)/2)+1):nrow(dff),] # make testing set
  
  # estimating the nuisance parameters
  PI_P_X_mu0 = glm(Y~X1+X2+X3+X4, data = dff.train[which(dff.train$A==0),], family = "poisson")
  PI_P_X_mu1 = glm(Y~X1+X2+X3+X4, data = dff.train[which(dff.train$A==1),], family = "binomial")
  
  # getting fitted values
  PI_P_X_mu0_hat = predict(PI_P_X_mu0, newdata = dff.test, type="response")
  PI_P_X_mu1_hat = predict(PI_P_X_mu1, newdata = dff.test, type="response")
  
  # calculating gamma hat
  PI_P_X_gammahat = (PI_P_X_mu1_hat - PI_P_X_mu0_hat)/PI_P_X_mu1_hat
  PI_P_X_gammahat_mean = mean(sample(PI_P_X_gammahat, 100))
  
  # RMSE
  PI_P_X_RMSE = sqrt(  mean( (PI_P_X_gammahat_mean - meangamma)^2 )  ) #
  
  # Bias
  PI_P_X_bias = mean(abs(PI_P_X_gammahat_mean - meangamma))
  
  # function: boostrap to estimate confidence interval
  fun.PI_P_X_boot = function(xthedata){
    
    # data splitting
    index = sample(1:nrow(xthedata)) # randomize indices
    xthedata.train = xthedata[index,][1:(nrow(xthedata)/2),] # make training set
    xthedata.test = xthedata[index,][((nrow(xthedata)/2)+1):nrow(xthedata),] # make testing set
    
    # estimating the nuisance parameters
    PI_P_X_mu0 = glm(Y~X1+X2+X3+X4, data = xthedata.train[which(xthedata.train$A==0),], family = "poisson")
    PI_P_X_mu1 = glm(Y~X1+X2+X3+X4, data = xthedata.train[which(xthedata.train$A==1),], family = "binomial")
    
    # getting fitted values 
    PI_P_X_mu0_hat = predict(PI_P_X_mu0, newdata = xthedata.test, type="response")
    PI_P_X_mu1_hat = predict(PI_P_X_mu1, newdata = xthedata.test, type="response")
    
    # calculating gamma hat
    PI_P_X_gammahat = (PI_P_X_mu1_hat - PI_P_X_mu0_hat)/PI_P_X_mu1_hat
    PI_P_X_gammahat_mean = mean(sample(PI_P_X_gammahat, 100))
    
    return(PI_P_X_gammahat_mean)
  }
  
  # apply bootstrap
  bootreps=1000
  PI_P_X_bootvec = 0
  #stime = system.time({ # this is for testing the time. It's 3 times faster with dopar than without.
  PI_P_X_bootvec.1 = foreach(i=1:bootreps, .options.multicore=list(preschedule=TRUE)) %dopar% {
    
    # randomize indices
    index.b = sample(1:nrow(dff), replace=FALSE) 
    
    # select new sample of rows from dff
    newdff = dff[index.b,] 
    
    # calculating gammahat_star (star because it's from a bootstrap)
    PI_P_X_gammahat_star_mean = fun.PI_P_X_boot(newdff) 
    
    # store in a vector
    PI_P_X_bootvec[i] <- PI_P_X_gammahat_star_mean
  }
  #})
  #stime
  
  # standard error
  PI_P_X_gammahat_sd = sd(unlist(PI_P_X_bootvec.1), na.rm = TRUE)
  
  # confidence interval
  PI_P_X_ci = paste(PI_P_X_gammahat_mean - 2*PI_P_X_gammahat_sd / sqrt(samplesize), 
                    PI_P_X_gammahat_mean + 2*PI_P_X_gammahat_sd / sqrt(samplesize), 
                    sep=", ")
  
  return(c(PI_P_X_RMSE, PI_P_X_bias, PI_P_X_ci))
}


  # 2. Plug-in, parametric, transformed X's (misspecified model)
fun.PI_P_Xstar = function(dff){
  
  # data splitting
  index = sample(1:nrow(dff)) # randomize indices
  dff.train = dff[index,][1:(nrow(dff)/2),] # make training set
  dff.test = dff[index,][((nrow(dff)/2)+1):nrow(dff),] # make testing set
  
  # estimating the nuisance parameters
  PI_P_Xstar_mu0 = glm(Y~X1star+X2star+X3star+X4star, data = dff.train[which(dff.train$A==0),], family = "binomial")
  PI_P_Xstar_mu1 = glm(Y~X1star+X2star+X3star+X4star, data = dff.train[which(dff.train$A==1),], family = "binomial")
  
  # getting fitted values
  PI_P_Xstar_mu0_hat = predict(PI_P_Xstar_mu0, newdata = dff.test, type="response")
  PI_P_Xstar_mu1_hat = predict(PI_P_Xstar_mu1, newdata = dff.test, type="response")
  
  # calculating gamma hat
  PI_P_Xstar_gammahat  = (PI_P_Xstar_mu1_hat - PI_P_Xstar_mu0_hat)/PI_P_Xstar_mu1_hat
  PI_P_Xstar_gammahat_mean = mean(sample(PI_P_Xstar_gammahat, size = 100))
  
  # RMSE 
  PI_P_Xstar_RMSE = sqrt(  mean( (PI_P_Xstar_gammahat_mean - meangamma)^2 )  )
  
  # bias
  PI_P_Xstar_bias = mean(abs(PI_P_Xstar_gammahat_mean - meangamma))
  
  # function: boostrap to estimate confidence interval
  fun.PI_P_Xstar_boot = function(xthedata){
    
    # data splitting
    index = sample(1:nrow(xthedata)) # randomize indices
    xthedata.train = xthedata[index,][1:(nrow(xthedata)/2),] # make training set
    xthedata.test = xthedata[index,][((nrow(xthedata)/2)+1):nrow(xthedata),] # make testing set
    
    # estimating the nuisance parameters
    PI_P_Xstar_mu0_hat = glm(Y~X1star+X2star+X3star+X4star, data = xthedata.train[which(xthedata.train$A==0),], family = "binomial")
    PI_P_Xstar_mu1_hat = glm(Y~X1star+X2star+X3star+X4star, data = xthedata.train[which(xthedata.train$A==1),], family = "binomial")
    
    # getting fitted values 
    PI_P_Xstar_mu0_hat = predict(PI_P_Xstar_mu0, newdata = xthedata.test, type="response")
    PI_P_Xstar_mu1_hat = predict(PI_P_Xstar_mu1, newdata = xthedata.test, type="response")
    
    # calculating gamma hat
    PI_P_Xstar_gammahat = (PI_P_Xstar_mu1_hat - PI_P_Xstar_mu0_hat)/PI_P_Xstar_mu1_hat
    PI_P_Xstar_gammahat_mean = mean(sample(PI_P_Xstar_gammahat, 100))
    
    return(PI_P_Xstar_gammahat_mean)
  }
  
  # apply bootstrap
  bootreps=1000
  PI_P_Xstar_bootvec = 0
  PI_P_Xstar_bootvec.1 = foreach(i=1:bootreps, .options.multicore=list(preschedule=TRUE)) %dopar% {
    
    # randomize indices
    index.b = sample(1:nrow(dff), replace=FALSE) 
    
    # select new sample of rows from dff
    newdff = dff[index.b,] 
    
    # calculating gammahat_star (star because it's from a bootstrap)
    PI_P_Xstar_gammahat_star_mean = fun.PI_P_Xstar_boot(newdff) 
    
    # store in a vector
    PI_P_Xstar_bootvec[i] <- PI_P_Xstar_gammahat_star_mean
  }
  
  # Standard error
  PI_P_Xstar_gammahat_sd = sd(unlist(PI_P_Xstar_bootvec.1), na.rm = TRUE)
  
  # Confidence interval
  PI_P_Xstar_ci = paste(PI_P_Xstar_gammahat_mean - 2*PI_P_Xstar_gammahat_sd / sqrt(samplesize), 
                        PI_P_Xstar_gammahat_mean + 2*PI_P_Xstar_gammahat_sd / sqrt(samplesize), 
                        sep=", ")

  return(c(PI_P_Xstar_RMSE, PI_P_Xstar_bias, PI_P_Xstar_ci))
}
  

  # 3. Plug-in, nonparametric, untransformed X's estimator:
fun.PI_N_X = function(dff){

  # data splitting
  index = sample(1:nrow(dff)) # randomize indices
  dff.train = dff[index,][1:(nrow(dff)/2),] # make training set
  dff.test = dff[index,][((nrow(dff)/2)+1):nrow(dff),] # make testing set
  
  # estimating the nuisance parameters
  PI_N_X_mu0 = randomForest(Y~X1+X2+X3+X4, data = dff.train[which(dff.train$A==0),], type=regression) #do.trace = 100, 
  PI_N_X_mu1 = randomForest(Y~X1+X2+X3+X4, data = dff.train[which(dff.train$A==1),], type=regression)
  
  # getting fitted values
  PI_N_X_mu0_hat = predict(PI_N_X_mu0, newdata = dff.test)
  PI_N_X_mu1_hat = predict(PI_N_X_mu1, newdata = dff.test)
  
  # calculating gamma hat
  PI_N_X_gammahat = (PI_N_X_mu1_hat - PI_N_X_mu0_hat)/PI_N_X_mu1_hat
  PI_N_X_gammahat_mean = mean(sample(PI_N_X_gammahat, 100))
  
  # RMSE 
  PI_N_X_RMSE = sqrt(  mean( (PI_N_X_gammahat_mean - meangamma)^2 )  )
  
  # Bias
  PI_N_X_bias = mean(abs(PI_N_X_gammahat_mean - meangamma))
  
  # Confidence interval
  # # No valid confidence interval here.
  PI_N_X_ci = NA
  
  return(c(PI_N_X_RMSE, PI_N_X_bias, PI_N_X_ci))
}
  

  # 4. Plug-in, nonparametric, transformed X's estimator:
fun.PI_N_Xstar = function(dff){

  # data splitting
  index = sample(1:nrow(dff)) # randomize indices
  dff.train = dff[index,][1:(nrow(dff)/2),] # make training set
  dff.test = dff[index,][((nrow(dff)/2)+1):nrow(dff),] # make testing set
  
  # estimating the nuisance parameters
  PI_N_Xstar_mu0 = randomForest(Y~X1star+X2star+X3star+X4star, data = dff.train[which(dff.train$A==0),], type=regression) #do.trace = 100, 
  PI_N_Xstar_mu1 = randomForest(Y~X1star+X2star+X3star+X4star, data = dff.train[which(dff.train$A==1),], type=regression)
  
  # getting fitted values
  PI_N_Xstar_mu0_hat = predict(PI_N_Xstar_mu0, newdata = dff.test)
  PI_N_Xstar_mu1_hat = predict(PI_N_Xstar_mu1, newdata = dff.test)
  
  # calculating gamma hat
  PI_N_Xstar_gammahat = (PI_N_Xstar_mu1_hat - PI_N_Xstar_mu0_hat)/PI_N_Xstar_mu1_hat
  PI_N_Xstar_gammahat_mean = mean(sample(PI_N_Xstar_gammahat, 100))
  
  # RMSE 
  PI_N_Xstar_RMSE = sqrt(  mean( (PI_N_Xstar_gammahat_mean - meangamma)^2 )  )
  
  # Bias
  PI_N_Xstar_bias = mean(abs(PI_N_Xstar_gammahat_mean - meangamma))
  
  # Confidence interval
  # # No valid confidence interval here.
  PI_N_Xstar_ci = NA
  
  return(c(PI_N_Xstar_RMSE, PI_N_Xstar_bias, PI_N_Xstar_ci))
}
  

  # 5. Influence-function-based estimator, estimate nuisance parameters with plugin parametric untransformed Xs ----
fun.IF_P_X = function(dff){

  # data splitting
  index = sample(1:nrow(dff)) # randomize indices
  dff.train = dff[index,][1:(nrow(dff)/2),] # make training set
  dff.test = dff[index,][((nrow(dff)/2)+1):nrow(dff),] # make testing set
  
  # estimating the nuisance parameters
  IF_P_X_mu0 = glm(Y~X1+X2+X3+X4, data = dff.train[which(dff.train$A==0),], family = "poisson")
  IF_P_X_mu1 = glm(Y~X1+X2+X3+X4, data = dff.train[which(dff.train$A==1),], family = "binomial")
  IF_P_X_pi = glm(A~X1+X2+X3+X4, data = dff.train, family = "binomial")
  
  # getting fitted values 
  IF_P_X_mu0_hat = predict(IF_P_X_mu0, newdata = dff.test, type="response")
  IF_P_X_mu1_hat = predict(IF_P_X_mu1, newdata = dff.test, type="response")
  IF_P_X_pi_hat = predict(IF_P_X_pi, newdata = dff.test, type="response")
  
  # defining my pseudo-outcome
  IF_P_X_ystar = (1/IF_P_X_mu1_hat)*((IF_P_X_mu0_hat/IF_P_X_mu1_hat)*(1/IF_P_X_pi_hat)*A*(Y-IF_P_X_mu1_hat) - 
                                       (1/(1-IF_P_X_pi_hat))*(1-A)*(Y-IF_P_X_mu0_hat)) + (IF_P_X_mu1_hat-IF_P_X_mu0_hat)/IF_P_X_mu1_hat
  
  # fitting model
  IF_P_X_model = try(
    nls(IF_P_X_ystar ~ fun.expit(beta1*X1 + beta2*X2 + beta3*X3 + beta4*X4),
        start=list(beta1=0, beta2=0, beta3=0, beta4=0),
        lower=c(-2, -2, -2, -2), 
        upper=c(2, 2, 2, 2),
        algorithm="port",
        data=dff, 
        nls.control(maxiter = 500)),
    silent=TRUE)
  
  # getting predicted values
  IF_P_X_gammahat = try(fun.expit(predict(IF_P_X_model)), silent=TRUE)
  IF_P_X_gammahat_mean = mean(sample(IF_P_X_gammahat, size = 100))
  
  # RMSE
  IF_P_X_RMSE = try(sqrt(  mean( (IF_P_X_gammahat_mean - meangamma)^2 )  ), silent=TRUE)
  
  # show NA if I get an error message
  IF_P_X_RMSE = ifelse(class(IF_P_X_RMSE)=="numeric", IF_P_X_RMSE, NA)
  
  # bias
  IF_P_X_bias = try( mean(abs(IF_P_X_gammahat_mean - meangamma)) , silent=TRUE)
  
  # show NA if I get an error message
  IF_P_X_bias = ifelse(class(IF_P_X_bias)=="numeric" && IF_P_X_bias!=Inf, IF_P_X_bias, NA)
  
  # function: boostrap to estimate confidence interval
  fun.IF_P_X_boot = function(xthedata){
    
    # data splitting
    index = sample(1:nrow(xthedata)) # randomize indices
    xthedata.train = xthedata[index,][1:(nrow(xthedata)/2),] # make training set
    xthedata.test = xthedata[index,][((nrow(xthedata)/2)+1):nrow(xthedata),] # make testing set
    
    # estimating the nuisance parameters
    IF_P_X_mu0 = glm(Y~X1+X2+X3+X4, data = xthedata.train[which(xthedata.train$A==0),], family = "poisson")
    IF_P_X_mu1 = glm(Y~X1+X2+X3+X4, data = xthedata.train[which(xthedata.train$A==1),], family = "binomial")
    IF_P_X_pi = glm(A~X1+X2+X3+X4, data = xthedata.train, family = "binomial")
    
    # getting fitted values 
    IF_P_X_mu0_hat = predict(IF_P_X_mu0, newdata = xthedata.test, type="response")
    IF_P_X_mu1_hat = predict(IF_P_X_mu1, newdata = xthedata.test, type="response")
    IF_P_X_pi_hat = predict(IF_P_X_pi, newdata = xthedata.test, type="response")
    
    # defining my pseudo-outcome
    IF_P_X_ystar = (1/IF_P_X_mu1_hat)*((IF_P_X_mu0_hat/IF_P_X_mu1_hat)*(1/IF_P_X_pi_hat)*A*(Y-IF_P_X_mu1_hat) - 
                                         (1/(1-IF_P_X_pi_hat))*(1-A)*(Y-IF_P_X_mu0_hat)) + (IF_P_X_mu1_hat-IF_P_X_mu0_hat)/IF_P_X_mu1_hat
    
    # fitting model
    IF_P_X_model = try(
      nls(IF_P_X_ystar ~ fun.expit(beta1*X1 + beta2*X2 + beta3*X3 + beta4*X4),
          start=list(beta1=0, beta2=0, beta3=0, beta4=0),
          lower=c(-2, -2, -2, -2),
          upper=c(2, 2, 2, 2),
          algorithm="port",
          data=xthedata, 
          nls.control(maxiter = 500)), 
      silent=TRUE)
    
    # getting predicted values for gamma hat
    IF_P_X_gammahat = try(fun.expit(predict(IF_P_X_model)), silent=TRUE)
    IF_P_X_gammahat_mean = mean(sample(IF_P_X_gammahat, size = 100))
    
    return(IF_P_X_gammahat_mean)
  }
  
  
  # apply bootstrap
  bootreps=10
  IF_P_X_bootvec = 0
  IF_P_X_bootvec.1 = foreach(i=1:bootreps, .options.multicore=list(preschedule=TRUE)) %dopar% {
    
    # randomize indices
    index.b = sample(1:nrow(dff), replace=FALSE) 
    
    # select new sample of rows from dff
    newdff = dff[index.b,] 
    
    # calculating gammahat_star (star because it's from a bootstrap)
    IF_P_X_gammahat_star_mean = fun.IF_P_X_boot(newdff) 
    
    # store in a vector
    IF_P_X_bootvec[i] <- IF_P_X_gammahat_star_mean
  }
  
  # standard error
  IF_P_X_gammahat_sd = sd(unlist(IF_P_X_bootvec.1), na.rm = TRUE)
  
  # confidence interval
  IF_P_X_ci = paste(IF_P_X_gammahat_mean - 2*IF_P_X_gammahat_sd / sqrt(samplesize), 
                    IF_P_X_gammahat_mean + 2*IF_P_X_gammahat_sd / sqrt(samplesize), 
                    sep=", ")
  
  return(c(IF_P_X_RMSE, IF_P_X_bias, IF_P_X_ci))
}
  

  # 6. Influence-function-based estimator, estimate nuisance parameters with plugin nonparametric untransformed Xs ----
fun.IF_N_X = function(dff){  

    # data splitting
  index = sample(1:nrow(dff)) # randomize indices
  dff.train = dff[index,][1:(nrow(dff)/2),] # make training set
  dff.test = dff[index,][((nrow(dff)/2)+1):nrow(dff),] # make testing set
  
  # estimating the nuisance parameters
  IF_N_X_mu0 = randomForest(Y~X1+X2+X3+X4, data = dff.train[which(dff.train$A==0),], type=regression) #do.trace = 100, 
  IF_N_X_mu1 = randomForest(Y~X1+X2+X3+X4, data = dff.train[which(dff.train$A==1),], type=regression)
  IF_N_X_pi = randomForest(A~X1+X2+X3+X4, data = dff.train, type=regression)
  
  # getting fitted values
  IF_N_X_mu0_hat = predict(IF_N_X_mu0, newdata = dff.test)
  IF_N_X_mu1_hat = predict(IF_N_X_mu1, newdata = dff.test)
  IF_N_X_pi_hat = predict(IF_N_X_pi, newdata = dff.test)
  
  # defining my pseudo-outcome
  IF_N_X_ystar = (1/IF_N_X_mu1_hat)*((IF_N_X_mu0_hat/IF_N_X_mu1_hat)*(1/IF_N_X_pi_hat)*A*(Y-IF_N_X_mu1_hat) - 
                                       (1/(1-IF_N_X_pi_hat))*(1-A)*(Y-IF_N_X_mu0_hat)) + (IF_N_X_mu1_hat-IF_N_X_mu0_hat)/IF_N_X_mu1_hat
  
  # fitting model
  IF_N_X_model = try(
    nls(IF_N_X_ystar ~ fun.expit(beta1*X1 + beta2*X2 + beta3*X3 + beta4*X4),
        start=list(beta1=0, beta2=0, beta3=0, beta4=0),
        lower=c(-2, -2, -2, -2),
        upper=c(2, 2, 2, 2),
        algorithm="port",
        data=dff, 
        nls.control(maxiter = 500)), 
    silent=TRUE)
  
  # getting predicted values
  IF_N_X_gammahat = try(fun.expit(predict(IF_N_X_model)), silent=TRUE)
  IF_N_X_gammahat_mean = mean(sample(IF_N_X_gammahat, size = 100))
  
  # RMSE
  IF_N_X_RMSE = try(sqrt(  mean( (IF_N_X_gammahat_mean - meangamma)^2 )  ), silent=TRUE)
  
  # show NA if I get an error message
  IF_N_X_RMSE = ifelse(class(IF_N_X_RMSE)=="numeric", IF_N_X_RMSE, NA)
  
  # bias
  IF_N_X_bias = try( mean(abs(IF_N_X_gammahat_mean - meangamma)) , silent=TRUE)
  
  # show NA if I get an error message
  IF_N_X_bias = ifelse(class(IF_N_X_bias)=="numeric" && IF_N_X_bias!=Inf, IF_N_X_bias, NA)
  
  # function: boostrap to estimate confidence interval
  fun.IF_N_X_boot = function(xthedata){
    
    # data splitting
    index = sample(1:nrow(xthedata)) # randomize indices
    xthedata.train = xthedata[index,][1:(nrow(xthedata)/2),] # make training set
    xthedata.test = xthedata[index,][((nrow(xthedata)/2)+1):nrow(xthedata),] # make testing set
    
    # estimating the nuisance parameters
    IF_N_X_mu0 = randomForest(Y~X1+X2+X3+X4, data = xthedata.train[which(xthedata.train$A==0),], type=regression) #do.trace = 100, 
    IF_N_X_mu1 = randomForest(Y~X1+X2+X3+X4, data = xthedata.train[which(xthedata.train$A==1),], type=regression)
    IF_N_X_pi = randomForest(A~X1+X2+X3+X4, data = xthedata.train, type=regression)
    
    # getting fitted values
    IF_N_X_mu0_hat = predict(IF_N_X_mu0, newdata = xthedata.test)
    IF_N_X_mu1_hat = predict(IF_N_X_mu1, newdata = xthedata.test)
    IF_N_X_pi_hat = predict(IF_N_X_pi, newdata = xthedata.test)
    
    # defining my pseudo-outcome
    IF_N_X_ystar = (1/IF_N_X_mu1_hat)*((IF_N_X_mu0_hat/IF_N_X_mu1_hat)*(1/IF_N_X_pi_hat)*A*(Y-IF_N_X_mu1_hat) - 
                                         (1/(1-IF_N_X_pi_hat))*(1-A)*(Y-IF_N_X_mu0_hat)) + (IF_N_X_mu1_hat-IF_N_X_mu0_hat)/IF_N_X_mu1_hat
    
    # fitting model
    IF_N_X_model = try(
      nls(IF_N_X_ystar ~ fun.expit(beta1*X1 + beta2*X2 + beta3*X3 + beta4*X4),
          start=list(beta1=0, beta2=0, beta3=0, beta4=0),
          lower=c(-2, -2, -2, -2),
          upper=c(2, 2, 2, 2),
          algorithm="port",
          data=xthedata, 
          nls.control(maxiter = 500)), 
      silent=TRUE)
    
    # getting predicted values
    IF_N_X_gammahat = try(fun.expit(predict(IF_N_X_model)), silent=TRUE)
    IF_N_X_gammahat_mean = mean(sample(IF_N_X_gammahat, size = 100))
    
    return(IF_N_X_gammahat_mean)
  }
  
  # apply bootstrap
  bootreps=10
  IF_N_X_bootvec = 0
  IF_N_X_bootvec.1 = foreach(i=1:bootreps, .options.multicore=list(preschedule=TRUE)) %dopar% {
    
    # randomize indices
    index.b = sample(1:nrow(dff), replace=FALSE) 
    
    # select new sample of rows from dff
    newdff = dff[index.b,] 
    
    # calculating gammahat_star (star because it's from a bootstrap)
    IF_N_X_gammahat_star_mean = fun.IF_N_X_boot(newdff) 
    
    # store in a vector
    IF_N_X_bootvec[i] <- IF_N_X_gammahat_star_mean
  }
  
  # standard error
  IF_N_X_gammahat_sd = sd(unlist(IF_N_X_bootvec.1), na.rm = TRUE)
  
  # confidence interval
  IF_N_X_ci = paste(IF_N_X_gammahat_mean - 2*IF_N_X_gammahat_sd / sqrt(samplesize), 
                    IF_N_X_gammahat_mean + 2*IF_N_X_gammahat_sd / sqrt(samplesize), 
                    sep=", ")
  
 return(c(IF_N_X_bias, IF_N_X_RMSE, IF_N_X_ci)) 
}
  

  # 7. Influence-function-based estimator, estimate nuisance parameters with plugin parametric transformed Xs ----
fun.IF_P_Xstar = function(dff){

  # data splitting
  index = sample(1:nrow(dff)) # randomize indices
  dff.train = dff[index,][1:(nrow(dff)/2),] # make training set
  dff.test = dff[index,][((nrow(dff)/2)+1):nrow(dff),] # make testing set
  
  # estimating the nuisance parameters
  IF_P_Xstar_mu0 = glm(Y~X1+X2+X3+X4, data = dff.train[which(dff.train$A==0),], family = "poisson")
  IF_P_Xstar_mu1 = glm(Y~X1+X2+X3+X4, data = dff.train[which(dff.train$A==1),], family = "binomial")
  IF_P_Xstar_pi = glm(A~X1+X2+X3+X4, data = dff.train, family = "binomial")
  
  # getting fitted values
  IF_P_Xstar_mu0_hat = predict(IF_P_Xstar_mu0, newdata = dff.test, type="response")
  IF_P_Xstar_mu1_hat = predict(IF_P_Xstar_mu1, newdata = dff.test, type="response")
  IF_P_Xstar_pi_hat = predict(IF_P_Xstar_pi, newdata = dff.test, type="response")
  
  # defining my pseudo-outcome
  IF_P_Xstar_ystar = (1/IF_P_Xstar_mu1_hat)*((IF_P_Xstar_mu0_hat/IF_P_Xstar_mu1_hat)*(1/IF_P_Xstar_pi_hat)*A*(Y-IF_P_Xstar_mu1_hat) - 
                                               (1/(1-IF_P_Xstar_pi_hat))*(1-A)*(Y-IF_P_Xstar_mu0_hat)) + (IF_P_Xstar_mu1_hat-IF_P_Xstar_mu0_hat)/IF_P_Xstar_mu1_hat
  
  # fitting model
  IF_P_Xstar_model = try(
    nls(IF_P_Xstar_ystar ~ fun.expit(beta1*X1star + beta2*X2star + beta3*X3star + beta4*X4star),
        start=list(beta1=0, beta2=0, beta3=0, beta4=0),
        data=dff,
        nls.control(maxiter = 500)),
    silent=TRUE)
  
  # getting predicted values
  IF_P_Xstar_gammahat = try(fun.expit(predict(IF_P_Xstar_model)), silent=TRUE)
  IF_P_Xstar_gammahat_mean = mean(sample(IF_P_Xstar_gammahat, size = 100))
  
  # RMSE
  IF_P_Xstar_RMSE = try(sqrt(  mean( (IF_P_Xstar_gammahat_mean - meangamma)^2 )  ), silent=TRUE)
  
  # show NA if I get an error message
  IF_P_Xstar_RMSE = ifelse(class(IF_P_Xstar_RMSE)=="numeric" && IF_P_Xstar_RMSE!=Inf, IF_P_Xstar_RMSE, NA)
  
  # bias
  IF_P_Xstar_bias = try( mean(abs(IF_P_Xstar_gammahat_mean - meangamma)) , silent=TRUE)
  
  # show NA if I get an error message
  IF_P_Xstar_bias = ifelse(class(IF_P_Xstar_bias)=="numeric" && IF_P_Xstar_bias!=Inf, IF_P_Xstar_bias, NA)
  
  # function: boostrap to estimate confidence interval
  fun.IF_P_Xstar_boot = function(xthedata){
    
    # data splitting
    index = sample(1:nrow(xthedata)) # randomize indices
    xthedata.train = xthedata[index,][1:(nrow(xthedata)/2),] # make training set
    xthedata.test = xthedata[index,][((nrow(xthedata)/2)+1):nrow(xthedata),] # make testing set
    
    # estimating the nuisance parameters
    IF_P_Xstar_mu0 = glm(Y~X1+X2+X3+X4, data = xthedata.train[which(xthedata.train$A==0),], family = "poisson")
    IF_P_Xstar_mu1 = glm(Y~X1+X2+X3+X4, data = xthedata.train[which(xthedata.train$A==1),], family = "binomial")
    IF_P_Xstar_pi = glm(A~X1+X2+X3+X4, data = xthedata.train, family = "binomial")
    
    # getting fitted values 
    IF_P_Xstar_mu0_hat = predict(IF_P_Xstar_mu0, newdata = xthedata.test, type="response")
    IF_P_Xstar_mu1_hat = predict(IF_P_Xstar_mu1, newdata = xthedata.test, type="response")
    IF_P_Xstar_pi_hat = predict(IF_P_Xstar_pi, newdata = xthedata.test, type="response")
    
    # defining my pseudo-outcome
    IF_P_Xstar_ystar = (1/IF_P_Xstar_mu1_hat)*((IF_P_Xstar_mu0_hat/IF_P_Xstar_mu1_hat)*(1/IF_P_Xstar_pi_hat)*A*(Y-IF_P_Xstar_mu1_hat) - 
                                                 (1/(1-IF_P_Xstar_pi_hat))*(1-A)*(Y-IF_P_Xstar_mu0_hat)) + (IF_P_Xstar_mu1_hat-IF_P_Xstar_mu0_hat)/IF_P_Xstar_mu1_hat
    
    # fitting model
    IF_P_Xstar_model = try(
      nls(IF_P_Xstar_ystar ~ fun.expit(beta1*X1star + beta2*X2star + beta3*X3star + beta4*X4star),
          start=list(beta1=0, beta2=0, beta3=0, beta4=0),
          lower=c(-2, -2, -2, -2),
          upper=c(2, 2, 2, 2),
          algorithm="port",
          data=xthedata, 
          nls.control(maxiter = 500)),
      silent=TRUE)
    
    # Getting predicted values
    IF_P_Xstar_gammahat = try(fun.expit(predict(IF_P_Xstar_model)), silent=TRUE)
    IF_P_Xstar_gammahat_mean = mean(sample(IF_P_Xstar_gammahat, size = 100))
    
    return(IF_P_Xstar_gammahat_mean)
  }
  
  # apply bootstrap
  bootreps=10
  IF_P_Xstar_bootvec = 0
  IF_P_Xstar_bootvec.1 = foreach(i=1:bootreps, .options.multicore=list(preschedule=TRUE)) %dopar% {
    
    # randomize indices
    index.b = sample(1:nrow(dff), replace=FALSE) 
    
    # select new sample of rows from dff
    newdff = dff[index.b,] 
    
    # calculating gammahat_star (star because it's from a bootstrap)
    IF_P_Xstar_gammahat_star_mean = fun.IF_P_Xstar_boot(newdff) 
    
    # store in a vector
    IF_P_Xstar_bootvec[i] <- IF_P_Xstar_gammahat_star_mean
  }
  
  # Standard error
  IF_P_Xstar_gammahat_sd = sd(unlist(IF_P_Xstar_bootvec.1), na.rm = TRUE)
  
  # Confidence interval
  IF_P_Xstar_ci = paste(IF_P_Xstar_gammahat_mean - 2*IF_P_Xstar_gammahat_sd / sqrt(samplesize), 
                        IF_P_Xstar_gammahat_mean + 2*IF_P_Xstar_gammahat_sd / sqrt(samplesize), 
                        sep=", ")
  
  return(c(IF_P_Xstar_bias, IF_P_Xstar_RMSE, IF_P_Xstar_ci))
 
}
  

  # 8. Influence-function-based estimator, estimate nuisance parameters with plugin nonparametric transformed Xs ----
fun.IF_N_Xstar = function(dff){

  # data splitting
  index = sample(1:nrow(dff)) # randomize indices
  dff.train = dff[index,][1:(nrow(dff)/2),] # make training set
  dff.test = dff[index,][((nrow(dff)/2)+1):nrow(dff),] # make testing set
  
  # estimating the nuisance parameters
  IF_N_Xstar_mu0 = randomForest(Y~X1star+X2star+X3star+X4star, data = dff.train[which(dff.train$A==0),], type=regression) #do.trace = 100, 
  IF_N_Xstar_mu1 = randomForest(Y~X1star+X2star+X3star+X4star, data = dff.train[which(dff.train$A==1),], type=regression)
  IF_N_Xstar_pi = randomForest(A~X1star+X2star+X3star+X4star, data = dff.train, type=regression)
  
  # getting fitted values
  IF_N_Xstar_mu0_hat = predict(IF_N_Xstar_mu0, newdata = dff.test)
  IF_N_Xstar_mu1_hat = predict(IF_N_Xstar_mu1, newdata = dff.test)
  IF_N_Xstar_pi_hat = predict(IF_N_Xstar_pi, newdata = dff.test)
  
  # defining my pseudo-outcome
  IF_N_Xstar_ystar = (1/IF_N_Xstar_mu1_hat)*((IF_N_Xstar_mu0_hat/IF_N_Xstar_mu1_hat)*(1/IF_N_Xstar_pi_hat)*A*(Y-IF_N_Xstar_mu1_hat) - 
                                               (1/(1-IF_N_Xstar_pi_hat))*(1-A)*(Y-IF_N_Xstar_mu0_hat)) + (IF_N_Xstar_mu1_hat-IF_N_Xstar_mu0_hat)/IF_N_Xstar_mu1_hat
  
  # fitting model
  IF_N_Xstar_model = try(
    nls(IF_N_Xstar_ystar ~ fun.expit(beta1*X1star + beta2*X2star + beta3*X3star + beta4*X4star),
        start=list(beta1=0, beta2=0, beta3=0, beta4=0),
        lower=c(-2, -2, -2, -2),
        upper=c(2, 2, 2, 2),
        algorithm="port",
        data=dff, 
        nls.control(maxiter = 500)),
    silent=TRUE)
  
  # getting predicted values
  IF_N_Xstar_gammahat = try(fun.expit(predict(IF_N_Xstar_model)), silent=TRUE)
  IF_N_Xstar_gammahat_mean = mean(sample(IF_N_Xstar_gammahat, size = 100))
  
  # RMSE
  IF_N_Xstar_RMSE = try(sqrt(  mean( (IF_N_Xstar_gammahat_mean - meangamma)^2 )  ), silent=TRUE)
  
  # show NA if I get an error message
  IF_N_Xstar_RMSE = ifelse(class(IF_N_Xstar_RMSE)=="numeric" && IF_N_Xstar_RMSE!=Inf, IF_N_Xstar_RMSE, NA)
  
  # bias
  IF_N_Xstar_bias = try( mean(abs(IF_N_Xstar_gammahat_mean - meangamma)) , silent=TRUE)
  
  # show NA if I get an error message
  IF_N_Xstar_bias = ifelse(class(IF_N_Xstar_bias)=="numeric" && IF_N_Xstar_bias!=Inf, IF_N_Xstar_bias, NA)
  
  # function: boostrap to estimate confidence interval
  fun.IF_N_Xstar_boot = function(xthedata){
    
    # data splitting
    index = sample(1:nrow(xthedata)) # randomize indices
    xthedata.train = xthedata[index,][1:(nrow(xthedata)/2),] # make training set
    xthedata.test = xthedata[index,][((nrow(xthedata)/2)+1):nrow(xthedata),] # make testing set
    
    # estimating the nuisance parameters
    IF_N_Xstar_mu0 = randomForest(Y~X1star+X2star+X3star+X4star, data = xthedata.train[which(xthedata.train$A==0),], type=regression) #do.trace = 100, 
    IF_N_Xstar_mu1 = randomForest(Y~X1star+X2star+X3star+X4star, data = xthedata.train[which(xthedata.train$A==1),], type=regression)
    IF_N_Xstar_pi = randomForest(A~X1star+X2star+X3star+X4star, data = xthedata.train, type=regression)
    
    # getting fitted values
    IF_N_Xstar_mu0_hat = predict(IF_N_Xstar_mu0, newdata = xthedata.test)
    IF_N_Xstar_mu1_hat = predict(IF_N_Xstar_mu1, newdata = xthedata.test)
    IF_N_Xstar_pi_hat = predict(IF_N_Xstar_pi, newdata = xthedata.test)
    
    # Defining my pseudo-outcome
    IF_N_Xstar_ystar = (1/IF_N_Xstar_mu1_hat)*((IF_N_Xstar_mu0_hat/IF_N_Xstar_mu1_hat)*(1/IF_N_Xstar_pi_hat)*A*(Y-IF_N_Xstar_mu1_hat) - 
                                                 (1/(1-IF_N_Xstar_pi_hat))*(1-A)*(Y-IF_N_Xstar_mu0_hat)) + (IF_N_Xstar_mu1_hat-IF_N_Xstar_mu0_hat)/IF_N_Xstar_mu1_hat
    
    # Fitting model
    IF_N_Xstar_model = try(
      nls(IF_N_Xstar_ystar ~ fun.expit(beta1*X1star + beta2*X2star + beta3*X3star + beta4*X4star),
          start=list(beta1=0, beta2=0, beta3=0, beta4=0),
          data=xthedata, 
          nls.control(maxiter = 500)),
      silent=TRUE)
    
    # Getting predicted values
    IF_N_Xstar_gammahat = try(fun.expit(predict(IF_N_Xstar_model)), silent=TRUE)
    IF_N_Xstar_gammahat_mean = mean(sample(IF_N_Xstar_gammahat, size = 100))
    
    return(IF_N_Xstar_gammahat_mean)
  }
  
  # apply bootstrap
  bootreps=10
  IF_N_Xstar_bootvec = 0
  IF_N_Xstar_bootvec.1 = foreach(i=1:bootreps, .options.multicore=list(preschedule=TRUE)) %dopar% {
    
    # randomize indices
    index.b = sample(1:nrow(dff), replace=FALSE) 
    
    # select new sample of rows from dff
    newdff = dff[index.b,] 
    
    # calculating gammahat_star (star because it's from a bootstrap)
    IF_N_Xstar_gammahat_star_mean = fun.IF_N_Xstar_boot(newdff) 
    
    # store in a vector
    IF_N_Xstar_bootvec[i] <- IF_N_Xstar_gammahat_star_mean
  }
  
  # Standard error
  IF_N_Xstar_gammahat_sd = sd(unlist(IF_N_Xstar_bootvec.1), na.rm = TRUE)
  
  # Confidence interval
  IF_N_Xstar_ci = paste(IF_N_Xstar_gammahat_mean - 2*IF_N_Xstar_gammahat_sd / sqrt(samplesize), 
                        IF_N_Xstar_gammahat_mean + 2*IF_N_Xstar_gammahat_sd / sqrt(samplesize), 
                        sep=", ")
  
  return(c(IF_N_Xstar_bias, IF_N_Xstar_RMSE, IF_N_Xstar_ci))
}



#----

# generate data 
dat = fun.generate.data(1000)

all.estimators = c(fun.PI_P_X, fun.PI_P_Xstar, fun.PI_N_X, fun.PI_N_Xstar,
                  fun.IF_P_X, fun.IF_N_X, fun.IF_P_Xstar, fun.IF_N_Xstar)

results = list()
results.list = list()
#stime = system.time({ # this is for testing the time. It's 3 times faster with dopar than without.
results = foreach(i=8, .options.multicore=list(preschedule=TRUE)) %dopar% {
    # store in a list
  results.list[[i]] = all.estimators[i][[1]](dat)
}
#})
#stime
results

fun.simulate(2000)
fun.simulate(200)
fun.simulate(200)

fun.simulate(500)
fun.simulate(500)
fun.simulate(500)
fun.simulate(500)

fun.simulate(1000)
fun.simulate(10000)



# # Repeat simulation for different sample sizes
reps = 10
samplesizes = c(rep(200, reps), rep(1000, reps), rep(10000, reps))
#samplesizes = c(rep(1000, 10), rep(2000, 10))


thelength = length(fun.simulate(1000))

arr <- array(dim=c(length(samplesizes),thelength+1))
colnames(arr) <- c("sample_sizes", 
                   "PI_P_X_RMSE", "PI_P_Xstar_RMSE", "PI_N_X_RMSE", "PI_N_Xstar_RMSE", 
                   "IF_P_X_RMSE", "IF_P_Xstar_RMSE", "IF_N_X_RMSE", "IF_N_Xstar_RMSE",
                   "PI_P_X_bias", "PI_P_Xstar_bias", "PI_N_X_bias", "PI_N_Xstar_bias", 
                   "IF_P_X_bias", "IF_P_Xstar_bias", "IF_N_X_bias", "IF_N_Xstar_bias",
                   "PI_P_X_ci", "PI_P_Xstar_ci", "PI_N_X_ci", "PI_N_Xstar_ci", 
                   "IF_P_X_ci", "IF_P_Xstar_ci", "IF_N_X_ci", "IF_N_Xstar_ci")


arr[1:length(samplesizes),1] = samplesizes
arr

for(s in 1:length(samplesizes)){
  fun.results = fun.simulate(samplesizes[s])
  for(j in 1:thelength){
    arr[s, j+1] = fun.results[j]
  }
}
arr
df.sim = as.data.frame(arr)


colnames(df.sim)[2] = "Plugin_param_cor_RMSE"
colnames(df.sim)[3] = "Plugin_param_mis_RMSE"
colnames(df.sim)[4] = "Plugin_nonp_cor_RMSE"
colnames(df.sim)[5] = "Plugin_nonp_mis_RMSE"
colnames(df.sim)[6] = "Influence_fun_param_cor_RMSE"
colnames(df.sim)[7] = "Influence_fun_param_mis_RMSE"
colnames(df.sim)[8] = "Influence_fun_nonp_cor_RMSE"
colnames(df.sim)[9] = "Influence_fun_nonp_mis_RMSE"
colnames(df.sim)[10] = "Plugin_param_cor_bias"
colnames(df.sim)[11] = "Plugin_param_mis_bias"
colnames(df.sim)[12] = "Plugin_nonp_cor_bias"
colnames(df.sim)[13] = "Plugin_nonp_mis_bias"
colnames(df.sim)[14] = "Influence_fun_param_cor_bias"
colnames(df.sim)[15] = "Influence_fun_param_mis_bias"
colnames(df.sim)[16] = "Influence_fun_nonp_cor_bias"
colnames(df.sim)[17] = "Influence_fun_nonp_mis_bias"
colnames(df.sim)[18] = "Plugin_param_cor_ci"
colnames(df.sim)[19] = "Plugin_param_mis_ci"
colnames(df.sim)[20] = "Plugin_nonp_cor_ci"
colnames(df.sim)[21] = "Plugin_nonp_mis_ci"
colnames(df.sim)[22] = "Influence_fun_param_cor_ci"
colnames(df.sim)[23] = "Influence_fun_param_mis_ci"
colnames(df.sim)[24] = "Influence_fun_nonp_cor_ci"
colnames(df.sim)[25] = "Influence_fun_nonp_mis_ci"



setwd(WD_simulation)
save(df.sim, file="df.sim.fin-2017-09-19_wcoverage.rda")
# load("df.sim.fin-2017-09-19_wcoverage.rda") # new one
# load("df.sim.fin-2017-05-19_rf.rda") # old one


# Estimating coverage
# - True value is mean(gamma)?? WHAT IS IT? It really should just be one value.


# Extract lower and upper bounds
fun.extract.lb = function(charvector){
  lb = substr(charvector, 1, regexpr(",", charvector)-1)
  return(lb)
}

fun.extract.ub = function(charvector){
  ub = substr(charvector, regexpr(",", charvector)+2, nchar(as.character(charvector[1])))
  return(ub)
}


# colnames(df.sim)[18] = "Plugin_param_cor_ci"
df.sim$Plugin_param_cor_ci_lb = as.numeric(fun.extract.lb(df.sim$Plugin_param_cor_ci))
df.sim$Plugin_param_cor_ci_ub = as.numeric(fun.extract.ub(df.sim$Plugin_param_cor_ci))

# colnames(df.sim)[19] = "Plugin_param_mis_ci"
df.sim$Plugin_param_mis_ci_lb = as.numeric(fun.extract.lb(df.sim$Plugin_param_mis_ci))
df.sim$Plugin_param_mis_ci_ub = as.numeric(fun.extract.ub(df.sim$Plugin_param_mis_ci))

# colnames(df.sim)[20] = "Plugin_nonp_cor_ci"
df.sim$Plugin_nonp_cor_ci_lb = as.numeric(fun.extract.lb(df.sim$Plugin_nonp_cor_ci))
df.sim$Plugin_nonp_cor_ci_ub = as.numeric(fun.extract.ub(df.sim$Plugin_nonp_cor_ci))

# colnames(df.sim)[21] = "Plugin_nonp_mis_ci"
df.sim$Plugin_nonp_mis_ci_lb = as.numeric(fun.extract.lb(df.sim$Plugin_nonp_mis_ci))
df.sim$Plugin_nonp_mis_ci_ub = as.numeric(fun.extract.ub(df.sim$Plugin_nonp_mis_ci))

# colnames(df.sim)[22] = "Influence_fun_param_cor_ci"
df.sim$Influence_fun_param_cor_ci_lb = as.numeric(fun.extract.lb(df.sim$Influence_fun_param_cor_ci))
df.sim$Influence_fun_param_cor_ci_ub = as.numeric(fun.extract.ub(df.sim$Influence_fun_param_cor_ci))

# colnames(df.sim)[23] = "Influence_fun_param_mis_ci"
df.sim$Influence_fun_param_mis_ci_lb = as.numeric(fun.extract.lb(df.sim$Influence_fun_param_mis_ci))
df.sim$Influence_fun_param_mis_ci_ub = as.numeric(fun.extract.ub(df.sim$Influence_fun_param_mis_ci))

# colnames(df.sim)[24] = "Influence_fun_nonp_cor_ci"
df.sim$Influence_fun_nonp_cor_ci_lb = as.numeric(fun.extract.lb(df.sim$Influence_fun_nonp_cor_ci))
df.sim$Influence_fun_nonp_cor_ci_ub = as.numeric(fun.extract.ub(df.sim$Influence_fun_nonp_cor_ci))

# colnames(df.sim)[25] = "Influence_fun_nonp_mis_ci"
df.sim$Influence_fun_nonp_mis_ci_lb = as.numeric(fun.extract.lb(df.sim$Influence_fun_nonp_mis_ci))
df.sim$Influence_fun_nonp_mis_ci_ub = as.numeric(fun.extract.ub(df.sim$Influence_fun_nonp_mis_ci))

# Keep only lower and upper bounds instead of ci variables
df.sim1 <- within(df.sim, rm("Plugin_param_cor_ci", "Plugin_param_mis_ci",
                             "Plugin_nonp_cor_ci", "Plugin_nonp_mis_ci", 
                             "Influence_fun_param_cor_ci", "Influence_fun_param_mis_ci", 
                             "Influence_fun_nonp_cor_ci", "Influence_fun_nonp_mis_ci"))


# These change if I decided to use different sample sizes
df.sim.200 <- df.sim1[which(df.sim1$sample_sizes==200),]
df.sim.1000 <- df.sim1[which(df.sim1$sample_sizes==1000),]
df.sim.10000 <- df.sim1[which(df.sim1$sample_sizes==10000),]


# True gamma?
meangamma = mean(gamma)
# FAKE GAMMA FOR NOW BECAUSE I CAN'T FIGURE IT OUT. I TRUST THE CI'S MORE THAN THE MEAN...
#meangamma = 0.627182 #(mean(df.sim.10000$Plugin_param_cor_lb)+mean(df.sim.10000$Plugin_param_cor_ub))*.5


# Calculating covarage
df.sims = c("df.sim.200", "df.sim.1000", "df.sim.10000")
thevec <- list()

for(i in 1:length(df.sims)){
  df.sim.n = get(df.sims[i])
  Plugin_param_cor_coverage = length(which(meangamma >= df.sim.n$Plugin_param_cor_ci_lb & meangamma <= df.sim.n$Plugin_param_cor_ci_ub))/nrow(df.sim.n)
  Plugin_param_mis_coverage   = length(which(meangamma >= df.sim.n$Plugin_param_mis_ci_lb & meangamma <= df.sim.n$Plugin_param_mis_ci_ub))/nrow(df.sim.n)
  Plugin_nonp_cor_coverage = length(which(meangamma >= df.sim.n$Plugin_nonp_cor_ci_lb & meangamma <= df.sim.n$Plugin_nonp_cor_ci_ub))/nrow(df.sim.n)
  Plugin_nonp_mis_coverage = length(which(meangamma >= df.sim.n$Plugin_nonp_mis_ci_lb & meangamma <= df.sim.n$Plugin_nonp_mis_ci_ub))/nrow(df.sim.n)
  Influence_fun_param_cor_coverage = length(which(meangamma >= df.sim.n$Influence_fun_param_cor_ci_lb & meangamma <= df.sim.n$Influence_fun_param_cor_ci_ub))/nrow(df.sim.n)
  Influence_fun_param_mis_coverage = length(which(meangamma >= df.sim.n$Influence_fun_param_mis_ci_lb & meangamma <= df.sim.n$Influence_fun_param_mis_ci_ub))/nrow(df.sim.n)
  Influence_fun_nonp_cor_coverage = length(which(meangamma >= df.sim.n$Influence_fun_nonp_cor_ci_lb & meangamma <= df.sim.n$Influence_fun_nonp_cor_ci_ub))/nrow(df.sim.n)
  Influence_fun_nonp_mis_coverage = length(which(meangamma >= df.sim.n$Influence_fun_nonp_mis_ci_lb & meangamma <= df.sim.n$Influence_fun_nonp_mis_ci_ub))/nrow(df.sim.n)
  
  thevec[[i]] <- c(Plugin_param_cor_coverage, Plugin_param_mis_coverage, 
                   Plugin_nonp_cor_coverage, Plugin_nonp_mis_coverage,
                   Influence_fun_param_cor_coverage, Influence_fun_param_mis_coverage, 
                   Influence_fun_nonp_cor_coverage, Influence_fun_nonp_mis_coverage)
  
}

thevec
length(unlist(thevec))


# Keep only coverage instead of lower and upper bounds
df.sim2 <- within(df.sim1, rm("Plugin_param_cor_ci_lb", "Plugin_param_cor_ci_ub", 
                              "Plugin_param_mis_ci_lb", "Plugin_param_mis_ci_ub",
                              "Plugin_nonp_cor_ci_lb", "Plugin_nonp_cor_ci_ub", 
                              "Plugin_nonp_mis_ci_lb", "Plugin_nonp_mis_ci_ub", 
                              "Influence_fun_param_cor_ci_lb", "Influence_fun_param_cor_ci_ub", 
                              "Influence_fun_param_mis_ci_lb", "Influence_fun_param_mis_ci_ub", 
                              "Influence_fun_nonp_cor_ci_lb", "Influence_fun_nonp_cor_ci_ub", 
                              "Influence_fun_nonp_mis_ci_lb", "Influence_fun_nonp_mis_ci_ub"))

df.sim2

df.sim2$Plugin_param_cor_coverage = c( rep(thevec[[1]][1], reps), rep(thevec[[2]][1], reps), rep(thevec[[3]][1], reps))
df.sim2$Plugin_param_mis_coverage = c( rep(thevec[[1]][2], reps), rep(thevec[[2]][2], reps), rep(thevec[[3]][2], reps))
df.sim2$Plugin_nonp_cor_coverage = NA
df.sim2$Plugin_nonp_mis_coverage = NA
df.sim2$Influence_fun_param_cor_coverage = c( rep(thevec[[1]][5], reps), rep(thevec[[2]][5], reps), rep(thevec[[3]][5], reps))
df.sim2$Influence_fun_param_mis_coverage = c( rep(thevec[[1]][6], reps), rep(thevec[[2]][6], reps), rep(thevec[[3]][6], reps))
df.sim2$Influence_fun_nonp_cor_coverage = c( rep(thevec[[1]][7], reps), rep(thevec[[2]][7], reps), rep(thevec[[3]][7], reps))
df.sim2$Influence_fun_nonp_mis_coverage = c( rep(thevec[[1]][8], reps), rep(thevec[[2]][8], reps), rep(thevec[[3]][8], reps))

length(df.sim2$sample_sizes)
length(df.sim2$Plugin_param_mis_coverage)

# ---------
# Newer draft with coverage
# ---------

# Barplot with error bars
# Reshape data frame from simulation, df.sim
aa = melt(data = df.sim2, id.vars = "sample_sizes", 
          variable.name = "Algorithm", value.name = "Value")
dat$Sample_Size = as.factor(dat$sample_sizes)

dat$RMSE = ifelse(str_detect(dat$Algorithm, "RMSE"), 1, 0)
dat$bias = ifelse(str_detect(dat$Algorithm, "bias"), 1, 0)
dat$coverage = ifelse(str_detect(dat$Algorithm, "coverage"), 1, 0)
dat$Value = as.numeric(dat$Value)

# Define function for mean and sd
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=FALSE),
      sd = sd(x[[col]], na.rm=FALSE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}



# Create new dataframe for barplot
dfsum = data_summary(data = dat, varname = "Value", groupnames = c("Sample_Size", "Algorithm"))
dfsum$Algorithm = as.factor(dfsum$Algorithm)
dfsum$theindex = 1:dim(dfsum)[1]
dfsum$Stat = ifelse(dfsum$theindex %in% grep("bias", dfsum$Algorithm), "Bias", 
                    ifelse(dfsum$theindex %in% grep("coverage", dfsum$Algorithm), "Coverage", "RMSE")
)
dfsum$Stat = as.factor(dfsum$Stat)
dfsum$Algorithm_type = ifelse(dfsum$theindex %in% grep("Plugin", dfsum$Algorithm), "Plugin", "Proposed")
dfsum$cormispnp = 0
dfsum$cormispnp = ifelse((str_detect(dfsum$Algorithm, "cor")& str_detect(dfsum$Algorithm, "param")), "Cor P", dfsum$cormispnp)
dfsum$cormispnp = ifelse((str_detect(dfsum$Algorithm, "mis")& str_detect(dfsum$Algorithm, "param")), "Mis P", dfsum$cormispnp)
dfsum$cormispnp = ifelse((str_detect(dfsum$Algorithm, "cor")& str_detect(dfsum$Algorithm, "nonp")), "Cor NP", dfsum$cormispnp)
dfsum$cormispnp = ifelse((str_detect(dfsum$Algorithm, "mis")& str_detect(dfsum$Algorithm, "nonp")), "Mis NP", dfsum$cormispnp)
dfsum$cormispnp = as.factor(dfsum$cormispnp)
dfsum$cormispnp = factor(dfsum$cormispnp,c("Cor P","Mis P","Cor NP", "Mis NP"))
dfsum$Sample_Sizes = paste("n=", dfsum$Sample_Size, sep="")
dfsum$Sample_Sizes = factor(dfsum$Sample_Sizes,c("n=200","n=1000", "n=10000")) # THIS CHANGES IF SAMPLE SIZES CHANGE
dfsum

head(dfsum)

bp = ggplot(dfsum, aes(x=cormispnp, y=Value, fill=Algorithm_type)) + 
  geom_bar(stat="identity", color="black", position="dodge", width=.5) +
  geom_errorbar(aes(ymin=Value-sd, ymax=Value+sd), width=.2, position = position_dodge(.5)) +
  ggtitle("Plug-in vs. proposed influence-function estimators, 10 iterations per sample size") + xlab("") + ylab("")


setwd(WD_figs)
pdf("20170920__Barplot.pdf", width=10, height=7)
bp + facet_grid(Stat ~ Sample_Sizes,scales = "free_y", switch = "y") + theme_bw()# + theme_minimal() 
dev.off()


















# ----------
# Older draft without coverage
# ----------



# Barplot with error bars
# Reshape data frame from simulatin, df.sim
aa = melt(data = df.sim, id.vars = "sample_sizes", 
          variable.name = "Algorithm", value.name = "RMSE") # Not changing the names right...
dat = rename(aa, c(sample_sizes = "Sample_Size", variable = "Algorithm", value = "RMSE"))
dat$Sample_Size = as.factor(dat$Sample_Size)
dat$Bias = c(rep(0,dim(dat)[1]/2), rep(1, dim(dat)[1]/2))



# Define function for mean and sd
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}




# Create new dataframe for barplot
dfsum = data_summary(data = dat, varname = "RMSE", groupnames = c("Sample_Size", "Algorithm"))
dfsum$Algorithm = as.factor(dfsum$Algorithm)
dfsum$theindex = 1:dim(dfsum)[1]
dfsum$bias = ifelse(dfsum$theindex %in% grep("bias", dfsum$Algorithm), "Bias", "RMSE")
dfsum$bias = as.factor(dfsum$bias)
#dfsum$bias = factor(dfsum$bias,c("RMSE","Bias"))
dfsum$Algorithm_type = ifelse(dfsum$theindex %in% grep("Plugin", dfsum$Algorithm), "Plugin", "Proposed")
dfsum$cormispnp = 0
dfsum$cormispnp = ifelse((str_detect(dfsum$Algorithm, "cor")& str_detect(dfsum$Algorithm, "param")), "Cor P", dfsum$cormispnp)
dfsum$cormispnp = ifelse((str_detect(dfsum$Algorithm, "mis")& str_detect(dfsum$Algorithm, "param")), "Mis P", dfsum$cormispnp)
dfsum$cormispnp = ifelse((str_detect(dfsum$Algorithm, "cor")& str_detect(dfsum$Algorithm, "nonp")), "Cor NP", dfsum$cormispnp)
dfsum$cormispnp = ifelse((str_detect(dfsum$Algorithm, "mis")& str_detect(dfsum$Algorithm, "nonp")), "Mis NP", dfsum$cormispnp)
dfsum$cormispnp = as.factor(dfsum$cormispnp)
dfsum$cormispnp = factor(dfsum$cormispnp,c("Cor P","Mis P","Cor NP", "Mis NP"))
#str_detect( str, "cor[\\s\\S]*param") # From
dfsum
dfsum$Sample_Sizes = paste("n=", dfsum$Sample_Size, sep="")
dfsum$Sample_Sizes = factor(dfsum$Sample_Sizes,c("n=200","n=1000", "n=10000")) # THIS CHANGES IF SAMPLE SIZES CHANGE

bp = ggplot(dfsum, aes(x=cormispnp, y=RMSE, fill=Algorithm_type)) + 
  geom_bar(stat="identity", color="black", position="dodge", width=.5) +
  geom_errorbar(aes(ymin=RMSE-sd, ymax=RMSE+sd), width=.2, position = position_dodge(.5)) +
  ggtitle("Plug-in vs. proposed influence-function estimators, 10 iterations per sample size") + xlab("") + ylab("")

bp

setwd(WD_figs)
pdf("20170519__Barplot_withbias.pdf", width=12, height=8)
bp + facet_grid(bias ~ Sample_Sizes,scales = "free_y", switch = "y") + theme_bw()# + theme_minimal() 
dev.off()






# Old barplot (without bias):

arr1 <- array(dim=c(length(samplesizes),thelength+1))
colnames(arr1) <- c("sample_sizes", "PI_P_X_RMSE", "PI_P_Xstar_RMSE", 
                    "PI_N_X_RMSE", "PI_N_Xstar_RMSE", "IF_P_X_RMSE", 
                    "IF_P_Xstar_RMSE", "IF_N_Xstar_RMSE")

arr1[1:length(samplesizes),1] = samplesizes
arr1


for(s in 1:length(samplesizes)){
  fun.results = fun.simulate(samplesizes[s])
  arr1[s, 2] = fun.results[1]
  arr1[s, 3] = fun.results[2]
  arr1[s, 4] = fun.results[3]
  arr1[s, 5] = fun.results[4]
  arr1[s, 6] = fun.results[5]
  arr1[s, 7] = fun.results[6]
  arr1[s, 8] = fun.results[7]
}
arr1
df.sim1 = as.data.frame(arr1)

colnames(df.sim1)[2] = "Plugin_param_cor"
colnames(df.sim1)[3] = "Plugin_param_mis"
colnames(df.sim1)[4] = "Plugin_nonp_cor"
colnames(df.sim1)[5] = "Plugin_nonp_mis"
colnames(df.sim1)[6] = "Influence_fun_param_cor"
colnames(df.sim1)[7] = "Influence_fun_param_mis"
colnames(df.sim1)[8] = "Influence_fun_nonp_cor"
colnames(df.sim1)[9] = "Influence_fun_nonp_mis"

setwd(WD_simulation)
save(df.sim1, file="df.sim.fin-2017-05-19_withoutbias.rda")

# # REMOVE OUTLIERS FOR Parametric_misspecified
# df.sim = df.sim[which(df.sim$Parametric_misspecified<0.5),]



# Alert when it's done running
beep(1) 



# Barplot with error bars (this is repeated from above)
# Reshape data frame from simulatin, df.sim
aa1 = melt(data = df.sim1, id.vars = "sample_sizes", 
           variable.name = "Algorithm", value.name = "RMSE") # Not changing the names right...
dat1 = rename(aa1, c(sample_sizes = "Sample_Size", variable = "Algorithm", value = "RMSE"))
dat1$Sample_Size = as.factor(dat1$Sample_Size)

# Define function for mean and sd
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}


# Create new dataframe for barplot
dfsum1 = data_summary(data = dat1, varname = "RMSE", groupnames = c("Sample_Size", "Algorithm"))
dfsum1$Algorithm = as.factor(dfsum1$Algorithm)


# Barplot with ggplot
bp1 = ggplot(dfsum, aes(x=Sample_Size, y=RMSE, fill=Algorithm)) + 
  geom_bar(stat="identity", color="black", position=position_dodge()) +
  geom_errorbar(aes(ymin=RMSE-sd, ymax=RMSE+sd), width=.2, position=position_dodge(.9)) +
  ggtitle("Plug-in vs. influence-function estimators, 500 iterations per sample size") + xlab("Sample size") + ylab("Root mean squared error")


setwd(WD_figs)
pdf("20170519_Barplot_withoutbias.pdf", width=8, height=3.3)
bp1 + scale_fill_brewer(palette="Paired")
dev.off()




# # Barplot with ggplot
# bp = ggplot(dfsum, aes(x=Algorithm, y=RMSE, fill=Algorithm)) + 
#   geom_bar(stat="identity", color="black", position=position_dodge()) +
#   geom_errorbar(aes(ymin=RMSE-sd, ymax=RMSE+sd), width=.2, position=position_dodge(.9)) +
#   ggtitle("Plug-in vs. influence-function estimators, 20 iterations per sample size") + xlab("Sample size") + ylab("Root mean squared error")

# setwd(WD_figs)
# pdf("20170517_long_Barplot.pdf", width=11, height=3.3)
# bp + scale_fill_brewer(palette="Paired")
# dev.off()


# Faceting with bias and sample size
bp + scale_fill_brewer(palette="Paired") + facet_grid(Bias ~ Sample_Size)





dfsum[,grep("PI", names(dfsum), value=TRUE)]







# Plotting the RMSE to compare estimators
mdf = melt(df.sim, id="sample_sizes")
mdf$Estimator = mdf$variable

setwd(WD_figs)
pdf(file="20170517_RMSE_comparison.pdf", width=10, height=5)
ggplot(na.omit(mdf), aes(x = sample_sizes, y=value, color=Estimator)) + geom_point() + 
  geom_smooth(method = "loess", size = 1.5) +
  labs(title = "Comparison of four estimators, 100 iterations per sample size", y="Root mean squared error", x="Sample size")
dev.off()



