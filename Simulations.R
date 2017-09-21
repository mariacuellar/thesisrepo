# # Thesis simulations
# # May 18, 2017
# # Making my own simulation for the parameter of interest: PC = gamma = 1-1/RR = 1-mu0/mu1.

# # define working directory
WD_figs = "/home/mcuellar"
WD_thesis = "/home/mcuellar/Thesis"
WD_simulation = "/home/mcuellar/Thesis"

# install.packages("np")
# install.packages("beepr")
# install.packages("randomForest")
# install.packages("ranger", dependencies=TRUE)
# install.packages("ggplot2")
# install.packages("tree")
# install.packages("SuperLearner")
# install.packages("polspline")
# install.packages("gbm")
# install.packages("glmnet")
# install.packages("earth")
# install.packages("wesanderson")
# install.packages("minpack.lm")
# install.packages("e1071")
# install.packages("reshape")
# install.packages("stringr")
# install.packages("ggthemes")
# install.packages('<R-package-file>')
# install.packages('boot')


# library(clusterPower) # lets me use expit function
# library(dplyr) # lets me use sample_n function
# library(hydroGOF) # lets me use rmse function
# library(np) # lets me use nonparametric models, like kernel (https://cran.r-project.org/web/packages/np/np.pdf)
# library(randomForest) # lets me do random forests
# library(ranger)# lets me do random forests
# library(e1071) # support vector machine
# library(neuralnet)
# library(reshape)
# library(ggplot2)
# library(tree)
# library(stats) # for nls
# library(SuperLearner)
# library(polspline)
# library(gbm)
# library(glmnet)
# library(earth)
# library(wesanderson)
# library(minpack.lm) # different type of nls
# library(e1071)
# library(reshape2)
# library(stringr)
# library('ggthemes')
# library(boot)

# Turn off warnings
#options(warn=-1)

# Turn on warnings
#options(warn=0)

#debug(fun.simulate)


# Function to simulate estimation of PC (gamma) as plugin, nonparametric, parametric IF, and nonparametric IF

samplesize = 2000 # test
fun.simulate = function(samplesize){
# set.seed(400) # seed for random number generator
# let us know how far it's gone
cat(paste(samplesize,"... "))

# true parameters
index = 1:samplesize
beta = 0.5
Y1 = rbinom(n = samplesize, size = 1, prob = beta)

X1 = rnorm(n = samplesize, mean = 0, sd = 1) # correct model
X2 = rnorm(n = samplesize, mean = 0, sd = 1)
X3 = rnorm(n = samplesize, mean = 0, sd = 1)
X4 = rnorm(n = samplesize, mean = 0, sd = 1)

X1star = exp(X1/2) # misspecified model
X2star = X2/(1 + exp(X1)) + 10
X3star = (X1*X3/25 + 0.6)^3
X4star = (X2 + X4 + 20)^2

pi = expit(-X1+0.5*X2-0.25*X3-0.1*X4)

A = rbinom(n = samplesize, size = 1, prob = expit(-X1+0.5*X2-0.25*X3-0.1*X4))

# Defining my parameter
beta = 0.5
mu0 = beta/(1+exp(-X1+0.5*X2-0.25*X3-0.1*X4))
mu1 = rep(beta, samplesize)
#gamma = 1 - mu0/mu1
gamma = expit(-X1+0.5*X2-0.25*X3-0.1*X4) # Don't we need to make it so that gamma is equal to 1-mu0/mu1?? It is.
meangamma = mean(sample(gamma, 100))

# generate data frame
df = as.data.frame(cbind(index, X1, X2, X3, X4, X1star, X2star, X3star, X4star, pi, gamma, A, Y1))
head(df)

# generate Y0 conditional on combinations of values of A and Y1
dfy11 = df[which(df$Y1==1),]
dfy11$Y0 = rbinom(n = nrow(dfy11) , size = 1, prob = (1-gamma)) # or is it location = expit(t(PSI)*X), scale = 0?

dfy10 = df[which(df$Y1==0),]
dfy10$Y0 = 0

# add Y0 to dataframe
df_wy0 = as.data.frame(rbind(dfy11, dfy10))

# apply consistency to get Y
df_wy0$Y = ifelse(df_wy0$A==1, df_wy0$Y1, df_wy0$Y0) 
Y = df_wy0$Y

# ordering data so it's as it was at the beginning
dff = df_wy0[order(df_wy0$index),] 
head(dff)




# now getting into the models ---

# Plug-in, parametric, untransformed X's (correctly specified model)
PI_P_X_mu0 = nls(Y~beta1/(1+exp(-X1+0.5*X2-0.25*X3-0.1*X4)), start=list(beta1=0.5), data = dff[which(dff$A==0),]) # take away the parameters and let it estimate them
PI_P_X_mu1 = glm(Y~X1+X2+X3+X4, data = dff[which(dff$A==1),], family = "binomial")

# getting fitted values (after inverse link function)
PI_P_X_mu0_hat = predict(PI_P_X_mu0, newdata = dff)
PI_P_X_mu1_hat = predict(PI_P_X_mu1, newdata = dff, type="response")

# calculating gamma hat
PI_P_X_gammahat = (PI_P_X_mu1_hat - PI_P_X_mu0_hat)/PI_P_X_mu1_hat

# RMSE
PI_P_X_RMSE = sqrt(  mean( (PI_P_X_gammahat - gamma)^2 )  ) #

# Bias
PI_P_X_bias = mean(abs(PI_P_X_gammahat - gamma))

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
for(i in 1:bootreps){
  index.b = sample(x = 1:nrow(dff), replace=TRUE) # randomize indices
  newdff = dff[index.b,] # select new sample of rows from dff
  PI_P_X_gammahat_stars = fun.PI_P_X_boot(newdff) # estimating gammahat for new dataset
  PI_P_X_gammahat_star = mean(sample(PI_P_X_gammahat_stars, size = 100)) # getting a sample of 100 gamma hat stars and taking their average
  PI_P_X_bootvec[i] <- PI_P_X_gammahat_star
  PI_P_X_bootvec
  }

# Standard error
PI_P_X_gammahat_sd = sd(PI_P_X_bootvec)

# Confidence interval
PI_P_X_ci_lb = mean(PI_P_X_gammahat) - 2*PI_P_X_gammahat_sd / sqrt(samplesize)
PI_P_X_ci_ub = mean(PI_P_X_gammahat) + 2*PI_P_X_gammahat_sd / sqrt(samplesize)
#PI_P_X_ci = c(PI_P_X_ci_lb, PI_P_X_ci_ub)
PI_P_X_ci = paste(PI_P_X_ci_lb, PI_P_X_ci_ub, sep=", ")






# Plug-in, parametric, transformed X's (misspecified model)
PI_P_Xstar_mu0 = glm(Y~X1star+X2star+X3star+X4star, data = dff[which(dff$A==0),], family = "binomial")
PI_P_Xstar_mu1 = glm(Y~X1star+X2star+X3star+X4star, data = dff[which(dff$A==1),], family = "binomial")

# getting fitted values (after inverse link function)
PI_P_Xstar_mu0_hat = predict(PI_P_Xstar_mu0, newdata = dff, type="response")
PI_P_Xstar_mu1_hat = predict(PI_P_Xstar_mu1, newdata = dff, type="response")

# calculating gamma hat
PI_P_Xstar_gammahat  = (PI_P_Xstar_mu1_hat - PI_P_Xstar_mu0_hat)/PI_P_Xstar_mu1_hat

# RMSE 
PI_P_Xstar_RMSE = sqrt(  mean( (PI_P_Xstar_gammahat - gamma)^2 )  ) # 0.2617908249
#sqrt(samplesize)*sqrt(  mean( (PI_P_Xstar_gammahat - gamma)^2 )  ) # 0.2617908249

# Bias
PI_P_Xstar_bias = mean(abs(PI_P_Xstar_gammahat - gamma))

# Coverage

# Coverage (Estimate coverage using the bootstrap)
# Function to estimate gamma
fun.PI_P_Xstar_boot = function(xthedata){
  # Plug-in, parametric, untransformed X's (correctly specified model)
  PI_P_Xstar_mu0 = glm(Y~X1star+X2star+X3star+X4star, data = xthedata[which(xthedata$A==0),], family = "binomial")
  PI_P_Xstar_mu1 = glm(Y~X1star+X2star+X3star+X4star, data = xthedata[which(xthedata$A==1),], family = "binomial")
  
  # getting fitted values (after inverse link function)
  PI_P_Xstar_mu0_hat = predict(PI_P_Xstar_mu0, newdata = xthedata, type="response")
  PI_P_Xstar_mu1_hat = predict(PI_P_Xstar_mu1, newdata = xthedata, type="response")
  
  # calculating gamma hat
  PI_P_Xstar_gammahat = (PI_P_Xstar_mu1_hat - PI_P_Xstar_mu0_hat)/PI_P_Xstar_mu1_hat
  return(PI_P_Xstar_gammahat)
}

# Bootstrap 
# bootreps defined above
PI_P_Xstar_bootvec <- 0
for(i in 1:bootreps){
  index.b = sample(x = 1:nrow(dff), replace=TRUE) # randomize indices
  newdff = dff[index.b,] # select new sample of rows from dff
  PI_P_Xstar_gammahat_stars = fun.PI_P_Xstar_boot(newdff) # estimating gammahat for new dataset
  PI_P_Xstar_gammahat_star = mean(sample(PI_P_Xstar_gammahat_stars, size = 100)) # getting a sample of 100 gamma hat stars and taking their average
  PI_P_Xstar_bootvec[i] <- PI_P_Xstar_gammahat_star
}

# Standard error
PI_P_Xstar_gammahat_sd = sd(PI_P_Xstar_bootvec)

# Confidence interval
PI_P_Xstar_ci_lb = mean(PI_P_Xstar_gammahat) - 2*PI_P_Xstar_gammahat_sd / sqrt(samplesize)
PI_P_Xstar_ci_ub = mean(PI_P_Xstar_gammahat) + 2*PI_P_Xstar_gammahat_sd / sqrt(samplesize)
#PI_P_Xstar_ci = c(PI_P_Xstar_ci_lb, PI_P_Xstar_ci_ub)
PI_P_Xstar_ci = paste(PI_P_Xstar_ci_lb, PI_P_Xstar_ci_ub, sep=", ")

# # Coverage
# PI_P_Xstar_coverage = length(which( dff$gamma >= PI_P_Xstar_ci_lb & dff$gamma <= PI_P_Xstar_ci_ub ))/samplesize






# Plug-in, nonparametric, untransformed X's estimator:
# RandomForest----
# Fitting model
PI_N_X_mu0 = randomForest(Y~X1+X2+X3+X4, data = dff[which(dff$A==0),], type=regression)
PI_N_X_mu1 = randomForest(Y~X1+X2+X3+X4, data = dff[which(dff$A==1),], type=regression)

# Getting fitted values (after inverse link function)
PI_N_X_mu0_hat = as.numeric(expit(predict(PI_N_X_mu0, newdata=dff))) # with random forest
PI_N_X_mu1_hat = as.numeric(expit(predict(PI_N_X_mu1, newdata=dff)))

# Gamma hat
PI_N_X_gammahat = (PI_N_X_mu1_hat - PI_N_X_mu0_hat)/PI_N_X_mu1_hat

# RMSE 
PI_N_X_RMSE = sqrt(  mean( (PI_N_X_gammahat - gamma)^2 )  )

# Bias
PI_N_X_bias = mean(abs(PI_N_X_gammahat - gamma))

# Confidence interval
# # No valid confidence interval here.
PI_N_X_ci = NA

# # Coverage
# # No valid confidence interval here.
# PI_N_X_coverage = NA



# # Support vector machine----
# # Fitting model:
# PI_N_X_mu0 = svm(Y~X1+X2+X3+X4, data = dff[which(dff$A==0),])
# PI_N_X_mu1 = svm(Y~X1+X2+X3+X4, data = dff[which(dff$A==1),])
# 
# # Getting predicted values
# PI_N_X_mu0_hat = expit(predict(PI_N_X_mu0, newdata = dff))
# PI_N_X_mu1_hat = expit(predict(PI_N_X_mu1, newdata = dff))
# 
# # Gamma hat
# PI_N_X_gammahat = (PI_N_X_mu1_hat - PI_N_X_mu0_hat)/PI_N_X_mu1_hat
# 
# # RMSE
# PI_N_X_RMSE = sqrt(  mean( (PI_N_X_gammahat - gamma)^2 )  )
# #sqrt(samplesize)*sqrt(  mean( (PI_N_X_gammahat - gamma)^2 )  )
# 
# # Bias
# PI_N_X_bias = mean(abs(PI_N_X_gammahat - gamma))





# Plug-in, nonparametric, transformed X's estimator:

# RandomForest----
# Fitting model
PI_N_Xstar_mu0 = randomForest(Y~X1+X2+X3+X4, data = dff[which(dff$A==0),], type=regression)
PI_N_Xstar_mu1 = randomForest(Y~X1+X2+X3+X4, data = dff[which(dff$A==1),], type=regression)

# Getting fitted values (after inverse link function)
PI_N_Xstar_mu0_hat = as.numeric(expit(predict(PI_N_Xstar_mu0, newdata=dff))) # with random forest
PI_N_Xstar_mu1_hat = as.numeric(expit(predict(PI_N_Xstar_mu1, newdata=dff)))

# Gamma hat
PI_N_Xstar_gammahat = (PI_N_Xstar_mu1_hat - PI_N_Xstar_mu0_hat)/PI_N_Xstar_mu1_hat

# RMSE 
PI_N_Xstar_RMSE = sqrt(  mean( (PI_N_Xstar_gammahat - gamma)^2 )  )

# Bias
PI_N_Xstar_bias = mean(abs(PI_N_Xstar_gammahat - gamma))

# Confidence interval
# No valid confidence interval here.
PI_N_Xstar_ci = NA

# # Coverage
# PI_N_Xstar_coverage = NA



# # Fitting model:
# PI_N_Xstar_mu0 = svm(Y~X1star+X2star+X3star+X4star, data = dff[which(dff$A==0),])
# PI_N_Xstar_mu1 = svm(Y~X1star+X2star+X3star+X4star, data = dff[which(dff$A==1),])
# 
# # Getting predicted values
# PI_N_Xstar_mu0_hat = expit(predict(PI_N_Xstar_mu0, newdata = dff))
# PI_N_Xstar_mu1_hat = expit(predict(PI_N_Xstar_mu1, newdata = dff))
# 
# # Gamma hat
# PI_N_Xstar_gammahat = (PI_N_Xstar_mu1_hat - PI_N_Xstar_mu0_hat)/PI_N_Xstar_mu1_hat
# 
# # RMSE
# PI_N_Xstar_RMSE = sqrt(  mean( (PI_N_Xstar_gammahat - gamma)^2 )  )
# #sqrt(samplesize)*sqrt(  mean( (PI_N_Xstar_gammahat - gamma)^2 )  )
# 
# # Bias
# PI_N_Xstar_bias = mean(abs(PI_N_Xstar_gammahat - gamma))
# 
# # Show NA if I get an error message
# PI_N_Xstar_bias = ifelse(class(PI_N_Xstar_bias)=="numeric" && PI_N_Xstar_bias!=Inf, PI_N_Xstar_bias, NA)




#  Influence-function-based estimator, estimate nuisance parameters with plugin parametric untransformed Xs ----
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
  IF_P_X_ystar ~ expit(beta1*X1 + beta2*X2 + beta3*X3 + beta4*X4),
  start=list(beta1=0, beta2=0, beta3=0, beta4=0),
  lower=c(-2, -2, -2, -2),
  upper=c(2, 2, 2, 2),
  algorithm="port",
  data=dff, nls.control(maxiter = 500))
  , silent=TRUE)


# Getting predicted values
IF_P_X_gammahat = try(predict(IF_P_X_model), silent=TRUE)

# RMSE
IF_P_X_RMSE = try(sqrt(  mean( (IF_P_X_gammahat - gamma)^2 )  ), silent=TRUE)
#try(sqrt(samplesize)*sqrt(  mean( (IF_P_X_gammahat - gamma)^2 )  ), silent=TRUE)

# Show NA if I get an error message
IF_P_X_RMSE = ifelse(class(IF_P_X_RMSE)=="numeric", IF_P_X_RMSE, NA)

# Bias
IF_P_X_bias = try( mean(abs(IF_P_X_gammahat - gamma)) , silent=TRUE)

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
  IF_P_X_mu0_hat = predict(IF_P_X_mu0, newdata = xthedata) # mu0
  IF_P_X_mu1_hat = predict(IF_P_X_mu1, newdata = xthedata, type="response") # mu1
  IF_P_X_pi_hat = predict(IF_P_X_pi, newdata = xthedata, type="response") #pi 
  
  # Defining my pseudo-outcome
  IF_P_X_ystar = (1/IF_P_X_mu1_hat)*((IF_P_X_mu0_hat/IF_P_X_mu1_hat)*(1/IF_P_X_pi_hat)*A*(Y-IF_P_X_mu1_hat) - 
                                       (1/(1-IF_P_X_pi_hat))*(1-A)*(Y-IF_P_X_mu0_hat)) + (IF_P_X_mu1_hat-IF_P_X_mu0_hat)/IF_P_X_mu1_hat
  
  # Fitting model
  IF_P_X_model = try(nls(
    IF_P_X_ystar ~ expit(beta1*X1 + beta2*X2 + beta3*X3 + beta4*X4),
    start=list(beta1=0, beta2=0, beta3=0, beta4=0),
    lower=c(-2, -2, -2, -2),
    upper=c(2, 2, 2, 2),
    algorithm="port",
    data=xthedata, nls.control(maxiter = 500))
    , silent=TRUE)
  
  # Getting predicted values for gamma hat
  IF_P_X_gammahat = try(predict(IF_P_X_model), silent=TRUE)
  
  return(IF_P_X_gammahat)
}

# Bootstrap
bootreps=1000
IF_P_X_bootvec = 0
for(i in 1:bootreps){
  index.b = sample(x = 1:nrow(dff), replace=TRUE) # randomize indices
  newdff = dff[index.b,] # select new sample of rows from dff
  IF_P_X_gammahat_stars = fun.IF_P_X_boot(newdff) # estimating gammahat for new dataset
  IF_P_X_gammahat_star = mean(sample(IF_P_X_gammahat_stars, size = 100)) # getting a sample of 100 gamma hat stars and taking their average
  IF_P_X_bootvec[i] <- IF_P_X_gammahat_star
}

# Standard error
IF_P_X_gammahat_sd = sd(IF_P_X_bootvec)

# Confidence interval
IF_P_X_ci_lb = mean(IF_P_X_gammahat) - 2*IF_P_X_gammahat_sd / sqrt(samplesize)
IF_P_X_ci_ub = mean(IF_P_X_gammahat) + 2*IF_P_X_gammahat_sd / sqrt(samplesize)
#IF_P_X_ci = c(IF_P_X_ci_lb, IF_P_X_ci_ub)
IF_P_X_ci = paste(IF_P_X_ci_lb, IF_P_X_ci_ub, sep=", ")






#  Influence-function-based estimator, estimate nuisance parameters with plugin nonparametric untransformed Xs ----
#  Defining the nuisance parameters
IF_N_X_mu0 = randomForest(Y~X1+X2+X3+X4, data = dff[which(dff$A==0),], type=regression) # same as PI_N_Xstar_mu0
IF_N_X_mu1 = randomForest(Y~X1+X2+X3+X4, data = dff[which(dff$A==1),], type=regression) # same as PI_N_Xstar_mu1
IF_N_X_pi = randomForest(A~X1+X2+X3+X4, data = dff, type=regression)  # pi # switch to pi if you want to cheat and do better

# Getting fitted values
IF_N_X_mu0_hat = as.numeric(expit(predict(IF_N_Xstar_mu0, newdata=dff))) # with random forest
IF_N_X_mu1_hat = as.numeric(expit(predict(IF_N_Xstar_mu1, newdata=dff)))
IF_N_X_pi_hat = as.numeric(expit(predict(IF_N_X_pi, newdata=dff)))

# Defining my pseudo-outcome
IF_N_X_ystar = (1/IF_N_X_mu1_hat)*((IF_N_X_mu0_hat/IF_N_X_mu1_hat)*(1/IF_N_X_pi_hat)*A*(Y-IF_N_X_mu1_hat) - 
                                 (1/(1-IF_N_X_pi_hat))*(1-A)*(Y-IF_N_X_mu0_hat)) + (IF_N_X_mu1_hat-IF_N_X_mu0_hat)/IF_N_X_mu1_hat

# Fitting model
IF_N_X_model = try(nls(
  IF_N_X_ystar ~ expit(beta1*X1 + beta2*X2 + beta3*X3 + beta4*X4),
  start=list(beta1=0, beta2=0, beta3=0, beta4=0),
  lower=c(-2, -2, -2, -2),
  upper=c(2, 2, 2, 2),
  algorithm="port",
  data=dff, nls.control(maxiter = 500))
  , silent=TRUE)

# Getting predicted values
IF_N_X_gammahat = try(predict(IF_N_X_model), silent=TRUE)

# RMSE
IF_N_X_RMSE = try(sqrt(  mean( (IF_N_X_gammahat - gamma)^2 )  ), silent=TRUE)
#try(sqrt(samplesize)*sqrt(  mean( (IF_N_X_gammahat - gamma)^2 )  ), silent=TRUE)

# Show NA if I get an error message
IF_N_X_RMSE = ifelse(class(IF_N_X_RMSE)=="numeric", IF_N_X_RMSE, NA)

# Bias
IF_N_X_bias = try( mean(abs(IF_N_X_gammahat - gamma)) , silent=TRUE)

# Show NA if I get an error message
IF_N_X_bias = ifelse(class(IF_N_X_bias)=="numeric" && IF_N_X_bias!=Inf, IF_N_X_bias, NA)

# Coverage (Estimate coverage using the bootstrap)
# Function to estimate gamma
fun.IF_N_X_boot = function(xthedata){
    #  Defining the nuisance parameters
    IF_N_X_mu0 = randomForest(Y~X1+X2+X3+X4, data = xthedata[which(xthedata$A==0),], type=regression) # same as PI_N_Xstar_mu0
    IF_N_X_mu1 = randomForest(Y~X1+X2+X3+X4, data = xthedata[which(xthedata$A==1),], type=regression) # same as PI_N_Xstar_mu1
    IF_N_X_pi = randomForest(A~X1+X2+X3+X4, data = xthedata, type=regression)  # pi # switch to pi if you want to cheat and do better
    
    # Getting fitted values
    IF_N_X_mu0_hat = as.numeric(expit(predict(PI_N_Xstar_mu0, newdata=xthedata))) # with random forest
    IF_N_X_mu1_hat = as.numeric(expit(predict(PI_N_Xstar_mu1, newdata=xthedata)))
    IF_N_X_pi_hat = as.numeric(expit(predict(IF_N_X_pi, newdata=xthedata)))
    
    # Defining my pseudo-outcome
    IF_N_X_ystar = (1/IF_N_X_mu1_hat)*((IF_N_X_mu0_hat/IF_N_X_mu1_hat)*(1/IF_N_X_pi_hat)*A*(Y-IF_N_X_mu1_hat) - 
                                         (1/(1-IF_N_X_pi_hat))*(1-A)*(Y-IF_N_X_mu0_hat)) + (IF_N_X_mu1_hat-IF_N_X_mu0_hat)/IF_N_X_mu1_hat
    
    # Fitting model
    IF_N_X_model = try(nls(
      IF_N_X_ystar ~ expit(beta1*X1 + beta2*X2 + beta3*X3 + beta4*X4),
      start=list(beta1=0, beta2=0, beta3=0, beta4=0),
      lower=c(-2, -2, -2, -2),
      upper=c(2, 2, 2, 2),
      algorithm="port",
      data=xthedata, nls.control(maxiter = 500))
      , silent=TRUE)
    
    # Getting predicted values
    IF_N_X_gammahat = try(predict(IF_N_X_model), silent=TRUE)
    
    return(IF_N_X_gammahat)
}

# Bootstrap
bootreps1=5
IF_N_X_bootvec = 0
for(i in 1:bootreps1){
    index.b = sample(x = 1:nrow(dff), replace=TRUE) # randomize indices
    newdff = dff[index.b,] # select new sample of rows from dff
    IF_N_X_gammahat_stars = fun.IF_N_X_boot(newdff) # estimating gammahat for new dataset
    IF_N_X_gammahat_star = mean(sample(IF_N_X_gammahat_stars, size = 100)) # getting a sample of 100 gamma hat stars and taking their average
    IF_N_X_bootvec[i] <- IF_N_X_gammahat_star
}

# Standard error
IF_N_X_gammahat_sd = sd(IF_N_X_bootvec)

# Confidence interval
IF_N_X_ci_lb = mean(IF_N_X_gammahat) - 2*IF_N_X_gammahat_sd / sqrt(samplesize)
IF_N_X_ci_ub = mean(IF_N_X_gammahat) + 2*IF_N_X_gammahat_sd / sqrt(samplesize)
#IF_N_X_ci = c(IF_N_X_ci_lb, IF_N_X_ci_ub)
IF_N_X_ci = paste(IF_N_X_ci_lb, IF_N_X_ci_ub, sep=", ")






#  Influence-function-based estimator, estimate nuisance parameters with plugin parametric transformed Xs ----
#  Defining the nuisance parameters
IF_P_Xstar_mu0 = nls(Y~beta1/(1+exp(-X1star+0.5*X2star-0.25*X3star-0.1*X4star)), start=list(beta1=0.5), data = dff[which(dff$A==0),]) # same as PI_P_X_mu0_hat
IF_P_Xstar_mu1 = glm(Y~X1star+X2star+X3star+X4star, data = dff[which(dff$A==1),], family = "binomial") # same as PI_P_X_mu1_hat
IF_P_Xstar_pi = glm(A~X1star+X2star+X3star+X4star, data = dff, family = "binomial") # same as PI_P_X_mu1_hat

# getting fitted values (after inverse link function)
IF_P_Xstar_mu0_hat = predict(IF_P_X_mu0, newdata = dff) # mu0 
IF_P_Xstar_mu1_hat = predict(IF_P_X_mu1, newdata = dff, type="response") # mu1
IF_P_Xstar_pi_hat = predict(IF_P_X_pi, newdata = dff, type="response") # pi

# Defining my pseudo-outcome
IF_P_Xstar_ystar = (1/IF_P_Xstar_mu1_hat)*((IF_P_Xstar_mu0_hat/IF_P_Xstar_mu1_hat)*(1/IF_P_Xstar_pi_hat)*A*(Y-IF_P_Xstar_mu1_hat) - 
                                     (1/(1-IF_P_Xstar_pi_hat))*(1-A)*(Y-IF_P_Xstar_mu0_hat)) + (IF_P_Xstar_mu1_hat-IF_P_Xstar_mu0_hat)/IF_P_Xstar_mu1_hat

# Fitting model
IF_P_Xstar_model = try(
  nls(IF_P_Xstar_ystar ~ expit(beta1*X1star + beta2*X2star + beta3*X3star + beta4*X4star),
      start=list(beta1=0, beta2=0, beta3=0, beta4=0),
      lower=c(-2, -2, -2, -2),
      upper=c(2, 2, 2, 2),
      algorithm="port",
      data=dff, nls.control(maxiter = 500))
    , silent=TRUE)

# Getting predicted values
IF_P_Xstar_gammahat = try(predict(IF_P_Xstar_model), silent=TRUE)

# RMSE
IF_P_Xstar_RMSE = try(sqrt(  mean( (IF_P_Xstar_gammahat - gamma)^2 )  ), silent=TRUE)
#try(sqrt(samplesize)*sqrt(  mean( (IF_P_Xstar_gammahat - gamma)^2 )  ), silent=TRUE)

# Show NA if I get an error message
IF_P_Xstar_RMSE = ifelse(class(IF_P_Xstar_RMSE)=="numeric" && IF_P_Xstar_RMSE!=Inf, IF_P_Xstar_RMSE, NA)

# Bias
IF_P_Xstar_bias = try( mean(abs(IF_P_Xstar_gammahat - gamma)) , silent=TRUE)

# Show NA if I get an error message
IF_P_Xstar_bias = ifelse(class(IF_P_Xstar_bias)=="numeric" && IF_P_Xstar_bias!=Inf, IF_P_Xstar_bias, NA)

# Coverage (Estimate coverage using the bootstrap)
# Function to estimate gamma
fun.IF_P_Xstar_boot = function(xthedata){
  #  Defining the nuisance parameters
  IF_P_Xstar_mu0 = nls(Y~beta1/(1+exp(-X1star+0.5*X2star-0.25*X3star-0.1*X4star)), start=list(beta1=0.5), data = xthedata[which(xthedata$A==0),]) # same as PI_P_X_mu0_hat
  IF_P_Xstar_mu1 = glm(Y~X1star+X2star+X3star+X4star, data = xthedata[which(xthedata$A==1),], family = "binomial") # same as PI_P_X_mu1_hat
  IF_P_Xstar_pi = glm(A~X1star+X2star+X3star+X4star, data = xthedata, family = "binomial") # same as PI_P_X_mu1_hat
  
  # getting fitted values (after inverse link function)
  IF_P_Xstar_mu0_hat = mu0 #predict(IF_P_Xstar_mu0, newdata = xthedata) # mu0 
  IF_P_Xstar_mu1_hat = mu1 #predict(IF_P_Xstar_mu1, newdata = xthedata, type="response") # mu1
  IF_P_Xstar_pi_hat = pi #predict(IF_P_Xstar_pi, newdata = xthedata, type="response") # pi
  
  # Defining my pseudo-outcome
  IF_P_Xstar_ystar = (1/IF_P_Xstar_mu1_hat)*((IF_P_Xstar_mu0_hat/IF_P_Xstar_mu1_hat)*(1/IF_P_Xstar_pi_hat)*A*(Y-IF_P_Xstar_mu1_hat) - 
                                               (1/(1-IF_P_Xstar_pi_hat))*(1-A)*(Y-IF_P_Xstar_mu0_hat)) + (IF_P_Xstar_mu1_hat-IF_P_Xstar_mu0_hat)/IF_P_Xstar_mu1_hat
  
  # Fitting model
  IF_P_Xstar_model = try(
    nls(IF_P_Xstar_ystar ~ expit(beta1*X1star + beta2*X2star + beta3*X3star + beta4*X4star),
        start=list(beta1=0, beta2=0, beta3=0, beta4=0),
        lower=c(-2, -2, -2, -2),
        upper=c(2, 2, 2, 2),
        algorithm="port",
        data=xthedata, nls.control(maxiter = 500))
    , silent=TRUE)
  
  # Getting predicted values
  IF_P_Xstar_gammahat = try(predict(IF_P_Xstar_model), silent=TRUE)
  
  return(IF_P_Xstar_gammahat)
}

# Bootstrap
bootreps=100
IF_P_Xstar_bootvec = 0
for(i in 1:bootreps){
  index.b = sample(x = 1:nrow(dff), replace=TRUE) # randomize indices
  newdff = dff[index.b,] # select new sample of rows from dff
  IF_P_Xstar_gammahat_stars = fun.IF_P_Xstar_boot(newdff) # estimating gammahat for new dataset
  IF_P_Xstar_gammahat_star = mean(sample(IF_P_Xstar_gammahat_stars, size = 100)) # getting a sample of 100 gamma hat stars and taking their average
  IF_P_Xstar_bootvec[i] <- IF_P_Xstar_gammahat_star
}

# Standard error
IF_P_Xstar_gammahat_sd = sd(IF_P_Xstar_bootvec)

# Confidence interval
IF_P_Xstar_ci_lb = mean(IF_P_Xstar_gammahat) - 2*IF_P_Xstar_gammahat_sd / sqrt(samplesize)
IF_P_Xstar_ci_ub = mean(IF_P_Xstar_gammahat) + 2*IF_P_Xstar_gammahat_sd / sqrt(samplesize)
#IF_P_Xstar_ci = c(IF_P_Xstar_ci_lb, IF_P_Xstar_ci_ub)
IF_P_Xstar_ci = paste(IF_P_Xstar_ci_lb, IF_P_Xstar_ci_ub, sep=", ")






#  Influence-function-based estimator, estimate nuisance parameters with plugin nonparametric transformed Xs ----
#  Defining the nuisance parameters
IF_N_Xstar_mu0 = randomForest(Y~X1star+X2star+X3star+X4star, data = dff[which(dff$A==0),], type=regression) # same as PI_N_Xstar_mu0
IF_N_Xstar_mu1 = randomForest(Y~X1star+X2star+X3star+X4star, data = dff[which(dff$A==1),], type=regression) # same as PI_N_Xstar_mu1
IF_N_Xstar_pi = randomForest(A~X1star+X2star+X3star+X4star, data = dff, type=regression)  # pi # switch to pi if you want to cheat and do better

# Getting fitted values
IF_N_Xstar_mu0_hat = as.numeric(expit(predict(IF_N_Xstar_mu0, newdata=dff))) # with random forest
IF_N_Xstar_mu1_hat = as.numeric(expit(predict(IF_N_Xstar_mu1, newdata=dff)))
IF_N_Xstar_pi_hat = as.numeric(expit(predict(IF_N_Xstar_pi, newdata=dff)))

# Defining my pseudo-outcome
IF_N_Xstar_ystar = (1/IF_N_Xstar_mu1_hat)*((IF_N_Xstar_mu0_hat/IF_N_Xstar_mu1_hat)*(1/IF_N_Xstar_pi_hat)*A*(Y-IF_N_Xstar_mu1_hat) - 
                                     (1/(1-IF_N_Xstar_pi_hat))*(1-A)*(Y-IF_N_Xstar_mu0_hat)) + (IF_N_Xstar_mu1_hat-IF_N_Xstar_mu0_hat)/IF_N_Xstar_mu1_hat

# Fitting model
IF_N_Xstar_model = try(nls(
  IF_N_Xstar_ystar ~ expit(beta1*X1star + beta2*X2star + beta3*X3star + beta4*X4star),
  start=list(beta1=0, beta2=0, beta3=0, beta4=0),
  lower=c(-2, -2, -2, -2),
  upper=c(2, 2, 2, 2),
  algorithm="port",
  data=dff, nls.control(maxiter = 500))
  , silent=TRUE)

# Getting predicted values
IF_N_Xstar_gammahat = try(predict(IF_N_Xstar_model), silent=TRUE)

# RMSE
IF_N_Xstar_RMSE = try(sqrt(  mean( (IF_N_Xstar_gammahat - gamma)^2 )  ), silent=TRUE)
#try(sqrt(samplesize)*sqrt(  mean( (IF_N_Xstar_gammahat - gamma)^2 )  ), silent=TRUE)

# Show NA if I get an error message
IF_N_Xstar_RMSE = ifelse(class(IF_N_Xstar_RMSE)=="numeric" && IF_N_Xstar_RMSE!=Inf, IF_N_Xstar_RMSE, NA)

# Bias
IF_N_Xstar_bias = try( mean(abs(IF_N_Xstar_gammahat - gamma)) , silent=TRUE)

# Show NA if I get an error message
IF_N_Xstar_bias = ifelse(class(IF_N_Xstar_bias)=="numeric" && IF_N_Xstar_bias!=Inf, IF_N_Xstar_bias, NA)

# Coverage
# Coverage (Estimate coverage using the bootstrap)
# Function to estimate gamma
fun.IF_N_Xstar_boot = function(xthedata){
  #  Defining the nuisance parameters
  IF_N_Xstar_mu0 = randomForest(Y~X1star+X2star+X3star+X4star, data = xthedata[which(xthedata$A==0),], type=regression) # same as PI_N_Xstar_mu0
  IF_N_Xstar_mu1 = randomForest(Y~X1star+X2star+X3star+X4star, data = xthedata[which(xthedata$A==1),], type=regression) # same as PI_N_Xstar_mu1
  IF_N_Xstar_pi = randomForest(A~X1star+X2star+X3star+X4star, data = xthedata, type=regression)  # pi # switch to pi if you want to cheat and do better
  
  # Getting fitted values
  IF_N_Xstar_mu0_hat = mu0 # as.numeric(expit(predict(IF_N_Xstar_mu0, newdata=xthedata))) # with random forest
  IF_N_Xstar_mu1_hat = mu1 # as.numeric(expit(predict(IF_N_Xstar_mu1, newdata=xthedata)))
  IF_N_Xstar_pi_hat = pi # as.numeric(expit(predict(IF_N_Xstar_pi, newdata=xthedata)))
  
  # Defining my pseudo-outcome
  IF_N_Xstar_ystar = (1/IF_N_Xstar_mu1_hat)*((IF_N_Xstar_mu0_hat/IF_N_Xstar_mu1_hat)*(1/IF_N_Xstar_pi_hat)*A*(Y-IF_N_Xstar_mu1_hat) - 
                                               (1/(1-IF_N_Xstar_pi_hat))*(1-A)*(Y-IF_N_Xstar_mu0_hat)) + (IF_N_Xstar_mu1_hat-IF_N_Xstar_mu0_hat)/IF_N_Xstar_mu1_hat
  
  # Fitting model
  IF_N_Xstar_model = try(
    nls(IF_N_Xstar_ystar ~ expit(beta1*X1star + beta2*X2star + beta3*X3star + beta4*X4star),
        start=list(beta1=0, beta2=0, beta3=0, beta4=0),
        lower=c(-2, -2, -2, -2),
        upper=c(2, 2, 2, 2),
        algorithm="port",
        data=xthedata, nls.control(maxiter = 500))
    , silent=TRUE)
  
  # Getting predicted values
  IF_N_Xstar_gammahat = try(predict(IF_N_Xstar_model), silent=TRUE)
  
  return(IF_N_Xstar_gammahat)
}

# Bootstrap
bootreps=100
IF_N_Xstar_bootvec = 0
for(i in 1:bootreps){
  index.b = sample(x = 1:nrow(dff), replace=TRUE) # randomize indices
  newdff = dff[index.b,] # select new sample of rows from dff
  IF_N_Xstar_gammahat_stars = fun.IF_N_Xstar_boot(newdff) # estimating gammahat for new dataset
  IF_N_Xstar_gammahat_star = mean(sample(IF_N_Xstar_gammahat_stars, size = 100)) # getting a sample of 100 gamma hat stars and taking their average
  IF_N_Xstar_bootvec[i] <- IF_N_Xstar_gammahat_star
}

# Standard error
IF_N_Xstar_gammahat_sd = sd(IF_N_Xstar_bootvec)

# Confidence interval
IF_N_Xstar_ci_lb = mean(IF_N_Xstar_gammahat) - 2*IF_N_Xstar_gammahat_sd / sqrt(samplesize)
IF_N_Xstar_ci_ub = mean(IF_N_Xstar_gammahat) + 2*IF_N_Xstar_gammahat_sd / sqrt(samplesize)
# IF_N_Xstar_ci = c(IF_N_Xstar_ci_lb, IF_N_Xstar_ci_ub)
IF_N_Xstar_ci = paste(IF_N_Xstar_ci_lb, IF_N_Xstar_ci_ub, sep=", ")




# Function will return this
toreturn = c(PI_P_X_RMSE, PI_P_Xstar_RMSE, PI_N_X_RMSE, PI_N_Xstar_RMSE, 
             IF_P_X_RMSE, IF_P_Xstar_RMSE, IF_N_X_RMSE, IF_N_Xstar_RMSE,
             PI_P_X_bias, PI_P_Xstar_bias, PI_N_X_bias, PI_N_Xstar_bias, 
             IF_P_X_bias, IF_P_Xstar_bias, IF_N_X_bias, IF_N_Xstar_bias,
             PI_P_X_ci, PI_P_Xstar_ci, PI_N_X_ci, PI_N_Xstar_ci, 
             IF_P_X_ci, IF_P_Xstar_ci, IF_N_X_ci, IF_N_Xstar_ci
             )

return(toreturn)
}



fun.simulate(200)
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






# FITTING OTHER NONPARAMETRIC PLUG IN ESTIMATORS



#   # GAM---- 
#   # Fitting model
#   cor_N_PI_mu0gam = gam(Y~X1+X2+X3+X4, family = binomial, data = dff[which(dff$A==0),])
#   cor_N_PI_mu1gam = gam(Y~X1+X2+X3+X4, family = binomial, data = dff[which(dff$A==1),])
#   
#   # Getting fitted values (after inverse link function)  
#   cor_N_PI_mu0_hatgam = expit(predict(cor_N_PI_mu0gam, newdata=dff))
#   cor_N_PI_mu1_hatgam = expit(predict(cor_N_PI_mu1gam, newdata=dff))
#   
#   # Gamma hat
#   cor_N_PI_gammahatgam = (cor_N_PI_mu1_hatgam - cor_N_PI_mu0_hatgam)/cor_N_PI_mu1_hatgam
#   
#   # RMSE
#   RMSE_cor_N_PIgam = sqrt(  mean( (cor_N_PI_gammahatgam - gamma)^2 )  )
#   RMSE_cor_N_PI = RMSE_cor_N_PIgam




#   # Ranger----
#   # Fitting model
#   cor_N_PI_mu0r = ranger(Y~X1+X2+X3+X4, data = dff[which(dff$A==0),])
#   cor_N_PI_mu1r = ranger(Y~X1+X2+X3+X4, data = dff[which(dff$A==1),])
#   
#   # Getting fitted values (after inverse link function)  
#   cor_N_PI_mu0_hatr = expit(predict(cor_N_PI_mu0r, data=dff)[1][[1]])
#   cor_N_PI_mu1_hatr = expit(predict(cor_N_PI_mu1r, data=dff)[1][[1]])
#   
#   # Gamma hat
#   cor_N_PI_gammahatr = (cor_N_PI_mu1_hatr - cor_N_PI_mu0_hatr)/cor_N_PI_mu1_hatr
#   
#   # RMSE
#   RMSE_cor_N_PIr = sqrt(  mean( (cor_N_PI_gammahatr - gamma)^2 )  )
# 
#   



# Misspecified model: NONparametric plug-in (P_PI) estimator:
# 
#   # Kernel----
#   # Fitting model, a bit faster than without tol, ftol, but still slow
#   # got bandwidths from running it once on a sample of the data
#   mis_N_PI_mu0k = npreg(Y~X1star+X2star+X3star+X4star, data = dff[which(dff$A==0),], tol=0.1, ftol=0.1, 
#                         bws=c(5229101.686, 590245.5023, 0.02962718249, 46.19266321))
#   
#   mis_N_PI_mu1k = npreg(Y~X1star+X2star+X3star+X4star, data = dff[which(dff$A==1),], tol=0.1, ftol=0.1, 
#                         bws=c(0.8807241131, 3242719.391, 0.0643948029, 85037300.36))
# 
#   # Getting fitted values
#   mis_N_PI_mu0_hatk = predict(mis_N_PI_mu0k, newdata = dff)#$mean # with faster kernel
#   mis_N_PI_mu1_hatk = predict(mis_N_PI_mu1k, newdata = dff)#$mean
#   
#   # Gamma hat
#   mis_N_PI_gammahatk = (mis_N_PI_mu1_hatk - mis_N_PI_mu0_hatk)/mis_N_PI_mu1_hatk
#   
#   # RMSE 
#   RMSE_mis_N_PIk = sqrt(  mean( (mis_N_PI_gammahatk - gamma)^2 )  ) # 0.2458764
#   RMSE_mis_N_PI = RMSE_mis_N_PIk


#   # GAM---- 
#   # Fitting model
#   mis_N_PI_mu0gam = gam(Y~X1star+X2star+X3star+X4star, family = binomial, data = dff[which(dff$A==0),])
#   mis_N_PI_mu1gam = gam(Y~X1star+X2star+X3star+X4star, family = binomial, data = dff[which(dff$A==1),])
#   
#   # Getting fitted values (after inverse link function)  
#   mis_N_PI_mu0_hatgam = expit(predict(mis_N_PI_mu0gam, newdata=dff))
#   mis_N_PI_mu1_hatgam = expit(predict(mis_N_PI_mu1gam, newdata=dff))
#   
#   # Gamma hat
#   mis_N_PI_gammahatgam = (mis_N_PI_mu1_hatgam - mis_N_PI_mu0_hatgam)/mis_N_PI_mu1_hatgam
#   
#   # RMSE
#   RMSE_mis_N_PIgam = sqrt(  mean( (mis_N_PI_gammahatgam - gamma)^2 )  )
#   RMSE_mis_N_PI = RMSE_mis_N_PIgam
#   



# # Ranger----
# # Fitting model
# mis_N_PI_mu0r = ranger(Y~X1star+X2star+X3star+X4star, data = dff[which(dff$A==0),])
# mis_N_PI_mu1r = ranger(Y~X1star+X2star+X3star+X4star, data = dff[which(dff$A==1),])
# 
# # Getting fitted values (after inverse link function)
# mis_N_PI_mu0_hatr = expit(predict(mis_N_PI_mu0r, data=dff)[1][[1]])
# mis_N_PI_mu1_hatr = expit(predict(mis_N_PI_mu1r, data=dff)[1][[1]])
# 
# # Gamma hat
# mis_N_PI_gammahatr = (mis_N_PI_mu1_hatr - mis_N_PI_mu0_hatr)/mis_N_PI_mu1_hatr
# 
# # RMSE
# RMSE_mis_N_PIr = sqrt(  mean( (mis_N_PI_gammahatr - gamma)^2 )  )





# # SuperLearner----
#   SL.ranger <- function (Y, X, newX, family, ...) { require("ranger")
#                                                   fit.rf <- ranger::ranger(Y ~ ., data=X); pred <- predict(fit.rf,data=newX)$predictions
#                                                   fit <- list(object = fit.rf); out <- list(pred = pred, fit = fit)
#                                                   class(out$fit) <- c("SL.ranger"); return(out) }
# 
#   sl.lib1 <- c("SL.earth","SL.gam","SL.gbm", "SL.glm","SL.glmnet","SL.mean")
# sl.lib2 <- c("SL.glm", "SL.randomForest", "SL.mean")
# #   sl.lib3 <- c("SL.glm", "SL.randomForest", "SL.gam", "SL.polymars", "SL.mean")
# 
# # Fitting model
# data0 = dff[which(dff$A==0),] # need to put data into superlearner as df
# data1 = dff[which(dff$A==1),]
# 
# mis_N_PI_mu0SL = SuperLearner(Y=data0$Y, X=as.data.frame(cbind(data0$X1star,data0$X2star,data0$X3star,data0$X4star)),
#                               SL.library = sl.lib2, family=binomial())
# mis_N_PI_mu1SL = SuperLearner(Y=data1$Y, X=as.data.frame(cbind(data1$X1star,data1$X2star,data1$X3star,data1$X4star)),
#                               SL.library = sl.lib2, family=binomial())
# 
# # Getting fitted values (after inverse link function?)
# mis_N_PI_mu0_hatSL = predict(mis_N_PI_mu0SL, newdata=as.data.frame(cbind(dff$X1star,dff$X2star,dff$X3star,dff$X4star)))$library.predict[,2]
# mis_N_PI_mu1_hatSL = predict(mis_N_PI_mu1SL, newdata=as.data.frame(cbind(dff$X1star,dff$X2star,dff$X3star,dff$X4star)))$library.predict[,2]
# 
# # Gamma hat
# mis_N_PI_gammahatSL = (mis_N_PI_mu1_hatSL - mis_N_PI_mu0_hatSL)/mis_N_PI_mu1_hatSL
# 
# # RMSE
# RMSE_mis_N_PISL = sqrt(  mean( (mis_N_PI_gammahatSL - gamma)^2 )  )
# RMSE_mis_N_PI = RMSE_mis_N_PISL



#   
#   # Tree ---- 
#   mis_N_PI_mu0tree = tree(factor(Y)~X1+X2+X3+X4, data = dff[which(dff$A==0),], )
#   mis_N_PI_mu1tree = tree(factor(Y)~X1+X2+X3+X4, data = dff[which(dff$A==1),])
#   
#   mis_N_PI_mu0_hattree = predict(mis_N_PI_mu0tree, newdata=dff)[,2]
#   mis_N_PI_mu1_hattree = predict(mis_N_PI_mu1tree, newdata=dff)[,2]
#   
#   mis_N_PI_gammahattree = (mis_N_PI_mu1_hattree - mis_N_PI_mu0_hattree)/(mis_N_PI_mu1_hattree)
#   
#   RMSE_mis_N_PItree = sqrt(  mean( (mis_N_PI_gammahattree - gamma)^2 )  ) # 0.288463088
#   RMSE_mis_N_PItree  


#   
# 
# RandomForest----
#   # Fitting model
#   mis_N_PI_mu0rf = randomForest(Y~X1+X2+X3+X4, data = dff[which(dff$A==0),], type=regression)
#   mis_N_PI_mu1rf = randomForest(Y~X1+X2+X3+X4, data = dff[which(dff$A==1),], type=regression)
#   
#   # Getting fitted values (after inverse link function)
#   mis_N_PI_mu0_hatrf = as.numeric(expit(predict(mis_N_PI_mu0rf, newdata=dff))) # with random forest
#   mis_N_PI_mu1_hatrf = as.numeric(expit(predict(mis_N_PI_mu1rf, newdata=dff)))
#   # 
#   # Gamma hat
#   mis_N_PI_gammahatrf = (mis_N_PI_mu1_hatrf - mis_N_PI_mu0_hatrf)/mis_N_PI_mu1_hatrf
#   
#   # RMSE 
#   RMSE_mis_N_PIrf = sqrt(  mean( (mis_N_PI_gammahatrf - gamma)^2 )  ) # 0.4826061
# #
#   
#   
#   
# # Support vector machine----
# # Fitting model:
# cor_N_PI_mu0svm = svm(Y~X1+X2+X3+X4, data = dff[which(dff$A==0),])
# cor_N_PI_mu1svm = svm(Y~X1+X2+X3+X4, data = dff[which(dff$A==1),])
# 
# # Getting predicted values
# cor_N_PI_mu0_hatsvm = expit(predict(cor_N_PI_mu0svm, newdata = dff))
# cor_N_PI_mu1_hatsvm = expit(predict(cor_N_PI_mu1svm, newdata = dff))
# 
# # Gamma hat
# cor_N_PI_gammahatsvm = (cor_N_PI_mu1_hatsvm - cor_N_PI_mu0_hatsvm)/cor_N_PI_mu1_hatsvm
# 
# # RMSE
# RMSE_cor_N_PIsvm = sqrt(  mean( (cor_N_PI_gammahatsvm - gamma)^2 )  ) # 0.3977806
# RMSE_cor_N_PI = RMSE_cor_N_PIsvm


# 
# 
# # Support vector machine----
# # Fitting model:
# mis_N_PI_mu0svm = svm(Y~X1star+X2star+X3star+X4star, data = dff[which(dff$A==0),])
# mis_N_PI_mu1svm = svm(Y~X1star+X2star+X3star+X4star, data = dff[which(dff$A==1),])
# 
# # Getting predicted values
# mis_N_PI_mu0_hatsvm = expit(predict(mis_N_PI_mu0svm, newdata = dff))
# mis_N_PI_mu1_hatsvm = expit(predict(mis_N_PI_mu1svm, newdata = dff))
# 
# # Gamma hat
# mis_N_PI_gammahatsvm = (mis_N_PI_mu1_hatsvm - mis_N_PI_mu0_hatsvm)/mis_N_PI_mu1_hatsvm
# 
# # RMSE
# RMSE_mis_N_PIsvm = sqrt(  mean( (mis_N_PI_gammahatsvm - gamma)^2 )  ) # 0.3977806
# RMSE_mis_N_PI = RMSE_mis_N_PIsvm
# 

# 
# # Ranger----
# # Fitting model
# cor_N_PI_mu0r = ranger(Y~X1+X2+X3+X4, data = dff[which(dff$A==0),])
# cor_N_PI_mu1r = ranger(Y~X1+X2+X3+X4, data = dff[which(dff$A==1),])
# 
# # Getting fitted values (after inverse link function)
# cor_N_PI_mu0_hatr = expit(predict(cor_N_PI_mu0r, data=dff)[1][[1]])
# cor_N_PI_mu1_hatr = expit(predict(cor_N_PI_mu1r, data=dff)[1][[1]])
# 
# # Gamma hat
# cor_N_PI_gammahatr = (cor_N_PI_mu1_hatr - cor_N_PI_mu0_hatr)/cor_N_PI_mu1_hatr
# 
# # RMSE
# RMSE_cor_N_PIr = sqrt(  mean( (cor_N_PI_gammahatr - gamma)^2 )  )
# RMSE_cor_N_PI = RMSE_cor_N_PIr
# 
# 
# 
# 
# # Misspecified model: NONparametric plug-in (P_PI) estimator:
# # Ranger----
# # Fitting model
# mis_N_PI_mu0r = ranger(Y~X1star+X2star+X3star+X4star, data = dff[which(dff$A==0),])
# mis_N_PI_mu1r = ranger(Y~X1star+X2star+X3star+X4star, data = dff[which(dff$A==1),])
# 
# # Getting fitted values (after inverse link function)
# mis_N_PI_mu0_hatr = expit(predict(mis_N_PI_mu0r, data=dff)[1][[1]])
# mis_N_PI_mu1_hatr = expit(predict(mis_N_PI_mu1r, data=dff)[1][[1]])
# 
# # Gamma hat
# mis_N_PI_gammahatr = (mis_N_PI_mu1_hatr - mis_N_PI_mu0_hatr)/mis_N_PI_mu1_hatr
# 
# # RMSE
# RMSE_mis_N_PIr = sqrt(  mean( (mis_N_PI_gammahatr - gamma)^2 )  )
# RMSE_mis_N_PI = RMSE_mis_N_PIr




# # Correctly specified model: NONparametric plug-in (P_PI) estimator:
# # RandomForest----
# # Fitting model
# cor_N_PI_mu0rf = randomForest(Y~X1+X2+X3+X4, data = dff[which(dff$A==0),], type=regression)
# cor_N_PI_mu1rf = randomForest(Y~X1+X2+X3+X4, data = dff[which(dff$A==1),], type=regression)
# 
# # Getting fitted values (after inverse link function)
# cor_N_PI_mu0_hatrf = as.numeric(expit(predict(cor_N_PI_mu0rf, newdata=dff))) # with random forest
# cor_N_PI_mu1_hatrf = as.numeric(expit(predict(cor_N_PI_mu1rf, newdata=dff)))
# 
# # Gamma hat
# cor_N_PI_gammahatrf = (cor_N_PI_mu1_hatrf - cor_N_PI_mu0_hatrf)/cor_N_PI_mu1_hatrf
# plot(cor_N_PI_gammahatrf, gamma)
# # RMSE
# RMSE_cor_N_PIrf = sqrt(samplesize)*sqrt(  mean( (cor_N_PI_gammahatrf - gamma)^2 )  ) # 0.4826061
# RMSE_cor_N_PI = RMSE_cor_N_PIrf
# 
# 
# 
# # Misspecified model: NONparametric plug-in (P_PI) estimator:
# # RandomForest----
# # Fitting model
# mis_N_PI_mu0rf = randomForest(Y~X1star+X2star+X3star+X4star, data = dff[which(dff$A==0),], type=regression)
# mis_N_PI_mu1rf = randomForest(Y~X1star+X2star+X3star+X4star, data = dff[which(dff$A==1),], type=regression)
# 
# # Getting fitted values (after inverse link function)
# mis_N_PI_mu0_hatrf = as.numeric(expit(predict(mis_N_PI_mu0rf, newdata=dff))) # with random forest
# mis_N_PI_mu1_hatrf = as.numeric(expit(predict(mis_N_PI_mu1rf, newdata=dff)))
# #
# # Gamma hat
# mis_N_PI_gammahatrf = (mis_N_PI_mu1_hatrf - mis_N_PI_mu0_hatrf)/mis_N_PI_mu1_hatrf
# 
# # RMSE
# RMSE_mis_N_PIrf = sqrt(samplesize)*sqrt(  mean( (mis_N_PI_gammahatrf - gamma)^2 )  ) # 0.4826061
# RMSE_mis_N_PI = RMSE_mis_N_PIrf





# # Kernel. Got bandwidths from running it once on a sample of the data
# PI_N_Xstar_mu0 = npreg(Y~X1star+X2star+X3star+X4star, data = dff[which(dff$A==0),], tol=0.1, ftol=0.1,
#                       bws=c(5229101.686, 590245.5023, 0.02962718249, 46.19266321))
# 
# PI_N_Xstar_mu1 = npreg(Y~X1star+X2star+X3star+X4star, data = dff[which(dff$A==1),], tol=0.1, ftol=0.1,
#                       bws=c(0.8807241131, 3242719.391, 0.0643948029, 85037300.36))
# 
# # Getting fitted values
# PI_N_Xstar_mu0_hat = predict(PI_N_Xstar_mu0, newdata = dff)#$mean # with faster kernel
# PI_N_Xstar_mu1_hat = predict(PI_N_Xstar_mu1, newdata = dff)#$mean
# 
# # Gamma hat
# PI_N_Xstar_gammahat = (PI_N_Xstar_mu1_hat - PI_N_Xstar_mu0_hat)/PI_N_Xstar_mu1_hat
# 
# # RMSE
# PI_N_Xstar_RMSE = sqrt(  mean( (PI_N_Xstar_gammahat - gamma)^2 )  ) # 0.2458764




# 
# # Kernel. Got bandwidths from running it once on a sample of the data
# PI_N_X_mu0 = npreg(Y~X1+X2+X3+X4, data = dff[which(dff$A==0),], tol=0.1, ftol=0.1,
#                       bws=c(10263069, 13280958, 0.3782153, 0.4474266))
# 
# PI_N_X_mu1 = npreg(Y~X1+X2+X3+X4, data = dff[which(dff$A==1),], tol=0.1, ftol=0.1,
#                       bws=c(2850532, 0.953212, 8485719, 0.7261236))
# 
# # Getting fitted values
# PI_N_X_mu0_hat = predict(PI_N_X_mu0, newdata = dff)
# PI_N_X_mu1_hat = predict(PI_N_X_mu1, newdata = dff)
# 
# # Gamma hat
# PI_N_X_gammahat = (PI_N_X_mu1_hat - PI_N_X_mu0_hat)/PI_N_X_mu1_hat
# 
# # RMSE
# PI_N_X_mu1_RMSE = sqrt(  mean( (PI_N_X_gammahat - gamma)^2 )  )
