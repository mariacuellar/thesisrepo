# # Thesis simulations
# # May 4, 2017
# # Making my own simulation for the parameter of interest: PC = gamma = 1-1/RR = 1-mu0/mu1.
# 
# # define working directory
# WD_figs = "/home/mcuellar"
# WD_thesis = "/home/mcuellar/Thesis"
# WD_simulation = "/home/mcuellar/Thesis"
# 
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


# Turn off warnings
#options(warn=-1)

# Turn on warnings
#options(warn=0)

#debug(fun.simulate)


# Function to simulate estimation of PC (gamma) as plugin, nonparametric, parametric IF, and nonparametric IF

samplesize = 200 # test
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
  gamma = expit(-X1+0.5*X2-0.25*X3-0.1*X4) # Don't we need to make it so that gamma is equal to 1-mu0/mu1??
  
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
  
  # Correctly specified model: parametric plug-in (P_PI) estimator:
  cor_P_PI_mu0 = nls(Y~beta1/(1+exp(-X1+0.5*X2-0.25*X3-0.1*X4)), start=list(beta1=0.5), data = dff[which(dff$A==0),]) # take away the parameters and let it estimate them
  cor_P_PI_mu1 = glm(Y~X1+X2+X3+X4, data = dff[which(dff$A==1),], family = "binomial")
  
  # getting fitted values (after inverse link function)
  cor_P_PI_mu0_hat = predict(cor_P_PI_mu0, newdata = dff)
  cor_P_PI_mu1_hat = predict.glm(cor_P_PI_mu1, newdata = dff, type="response")
  
  # calculating gamma hat
  cor_P_PI_gammahat = (cor_P_PI_mu1_hat - cor_P_PI_mu0_hat)/cor_P_PI_mu1_hat
  
  # RMSE 
  RMSE_cor_P_PI = sqrt(samplesize)*sqrt(  mean( (cor_P_PI_gammahat - gamma)^2 )  ) # 0.310089068
  
  
  
  
  # Misspecified model: parametric plug-in (P_PI) estimator:
  mis_P_PI_mu0 = glm(Y~X1star+X2star+X3star+X4star, data = dff[which(dff$A==0),], family = "binomial")
  mis_P_PI_mu1 = glm(Y~X1star+X2star+X3star+X4star, data = dff[which(dff$A==1),], family = "binomial")
  
  # getting fitted values (after inverse link function)
  mis_P_PI_mu0_hat = predict.glm(mis_P_PI_mu0, newdata = dff, type="response")
  mis_P_PI_mu1_hat = predict.glm(mis_P_PI_mu1, newdata = dff, type="response")
  
  # calculating gamma hat
  mis_P_PI_gammahat = (mis_P_PI_mu1_hat - mis_P_PI_mu0_hat)/mis_P_PI_mu1_hat
  
  # RMSE 
  RMSE_mis_P_PI = sqrt(samplesize)*sqrt(  mean( (mis_P_PI_gammahat - gamma)^2 )  ) # 0.2617908249
  
  
  
  
  # Correctly specified model: NONparametric plug-in (P_PI) estimator:
  # RandomForest----
  # Fitting model
  cor_N_PI_mu0rf = randomForest(Y~X1+X2+X3+X4, data = dff[which(dff$A==0),], type=regression)
  cor_N_PI_mu1rf = randomForest(Y~X1+X2+X3+X4, data = dff[which(dff$A==1),], type=regression)

  # Getting fitted values (after inverse link function)
  cor_N_PI_mu0_hatrf = as.numeric(expit(predict(cor_N_PI_mu0rf, newdata=dff))) # with random forest
  cor_N_PI_mu1_hatrf = as.numeric(expit(predict(cor_N_PI_mu1rf, newdata=dff)))
  
  # Gamma hat
  cor_N_PI_gammahatrf = (cor_N_PI_mu1_hatrf - cor_N_PI_mu0_hatrf)/cor_N_PI_mu1_hatrf

  # RMSE
  RMSE_cor_N_PIrf = sqrt(samplesize)*sqrt(  mean( (cor_N_PI_gammahatrf - gamma)^2 )  ) # 0.4826061
  RMSE_cor_N_PI = RMSE_cor_N_PIrf
  
  
  
  # Misspecified model: NONparametric plug-in (P_PI) estimator:
  # RandomForest----
  # Fitting model
  mis_N_PI_mu0rf = randomForest(Y~X1star+X2star+X3star+X4star, data = dff[which(dff$A==0),], type=regression)
  mis_N_PI_mu1rf = randomForest(Y~X1star+X2star+X3star+X4star, data = dff[which(dff$A==1),], type=regression)
  
  # Getting fitted values (after inverse link function)
  mis_N_PI_mu0_hatrf = as.numeric(expit(predict(mis_N_PI_mu0rf, newdata=dff))) # with random forest
  mis_N_PI_mu1_hatrf = as.numeric(expit(predict(mis_N_PI_mu1rf, newdata=dff)))
  #
  # Gamma hat
  mis_N_PI_gammahatrf = (mis_N_PI_mu1_hatrf - mis_N_PI_mu0_hatrf)/mis_N_PI_mu1_hatrf
  
  # RMSE
  RMSE_mis_N_PIrf = sqrt(samplesize)*sqrt(  mean( (mis_N_PI_gammahatrf - gamma)^2 )  ) # 0.4826061
  RMSE_mis_N_PI = RMSE_mis_N_PIrf

  
  
  
  #  Influence-function-based estimator, untransformed Xs ----
  #  Defining the nuisance parameters
  cor_IF_mu0 = cor_P_PI_mu0_hat
  cor_IF_mu1 = cor_P_PI_mu1_hat
  cor_IF_pi = pi #as.numeric(expit(predict(glm(A~X1+X2+X3+X4, data = dff, family = "binomial"))))
  
  # Defining my pseudo-outcome, which will be the "y" in my glm model, from my influence function
  cor_ystar = (1/cor_IF_mu1)*((cor_IF_mu0/cor_IF_mu1)*(1/cor_IF_pi)*A*(Y-cor_IF_mu1) - 
                                (1/(1-cor_IF_pi))*(1-A)*(Y-cor_IF_mu0)) + (cor_IF_mu1-cor_IF_mu0)/cor_IF_mu1
  
  # Fitting model
  # cropping percentiles from top and bottom
  dff$cor_ystar = cor_ystar
  quantile_lo_cor = quantile(cor_ystar, .005)[[1]]
  quantile_hi_cor = quantile(cor_ystar, .995)[[1]]

  dff00 = dff
  dff00$cor_ystar = ifelse(dff00$cor_ystar>quantile_lo_cor, dff00$cor_ystar, quantile_lo_cor)
  dff00$cor_ystar = ifelse(dff00$cor_ystar<quantile_hi_cor, dff00$cor_ystar, quantile_hi_cor)
  
  cor_model_IF = try(nls(
    cor_ystar ~ expit(beta1*X1 + beta2*X2 + beta3*X3 + beta4*X4),
    start=list(beta1=-1, beta2=0.5, beta3=-0.25, beta4=0.1),
    data=dff00, nls.control(maxiter = 500)
  ), silent=TRUE)

  # cor_model_IF = try(nls(
  #   cor_ystar ~ expit(beta1*X1 + beta2*X2 + beta3*X3 + beta4*X4),
  #   start=list(beta1=-1, beta2=0.5, beta3=-0.25, beta4=0.1),
  #   data=dff, nls.control(maxiter = 500)
  # ), silent=TRUE)
  
  # Getting predicted values
  cor_gammahat_IF = try(predict(cor_model_IF), silent=TRUE)
  
  # RMSE
  RMSE_cor_IF = try(sqrt(samplesize)*sqrt(  mean( (cor_gammahat_IF - gamma)^2 )  ), silent=TRUE)
  
  # Show NA if I get an error message
  RMSE_cor_IF = ifelse(class(RMSE_cor_IF)=="numeric", RMSE_cor_IF, NA)
  
  
  
  
  #  Influence-function-based estimator, transformed Xs ----
  #  Defining the nuisance parameters
  mis_IF_mu0 = mis_P_PI_mu0_hat
  mis_IF_mu1 = mis_P_PI_mu1_hat
  mis_IF_pi = pi #as.numeric(expit(predict(glm(A~X1star+X2star+X3star+X4star, data = dff, family = "binomial"))))

  # Defining my pseudo-outcome, which will be the "y" in my glm model, from my influence function
  mis_ystar = (1/mis_IF_mu1)*((mis_IF_mu0/mis_IF_mu1)*(1/mis_IF_pi)*A*(Y-mis_IF_mu1) -
                                (1/(1-mis_IF_pi))*(1-A)*(Y-mis_IF_mu0)) + (mis_IF_mu1-mis_IF_mu0)/mis_IF_mu1

  # Fitting model
  # cropping percentiles from top and bottom
  dff$mis_ystar = mis_ystar
  quantile_lo_mis = quantile(mis_ystar, .005)[[1]]
  quantile_hi_mis = quantile(mis_ystar, .995)[[1]]

  dff11 = dff
  dff11$mis_ystar = ifelse(dff11$mis_ystar<quantile_lo_mis, quantile_lo_mis, dff11$mis_ystar)
  dff11$mis_ystar = ifelse(dff11$cor_ystar>quantile_hi_mis, quantile_hi_mis, dff11$mis_ystar)

  mis_model_IF = try(nls(
    mis_ystar ~ expit(beta1*X1star + beta2*X2star + beta3*X3star + beta4*X4star),
    start=list(beta1=0, beta2=0, beta3=0, beta4=0),
    data=dff11, nls.control(maxiter = 500)
    ), silent=TRUE)

  # mis_model_IF = try(nls(
  #   mis_ystar ~ expit(beta1*X1star + beta2*X2star + beta3*X3star + beta4*X4star),
  #   start=list(beta1=0, beta2=0, beta3=0, beta4=0),
  #   data=dff, nls.control(maxiter = 500)
  # ), silent=TRUE)

  # Getting predicted values
  mis_gammahat_IF = try(predict(mis_model_IF), silent=TRUE)

  # RMSE
  RMSE_mis_IF = try(sqrt(samplesize)*sqrt(  mean( (mis_gammahat_IF - gamma)^2 )  ), silent=TRUE)
  
  # Show NA if I get an error message
  RMSE_mis_IF = ifelse(class(RMSE_mis_IF)=="numeric", RMSE_mis_IF, NA)

  
  
  # Function will return this
  toreturn = c(RMSE_cor_P_PI, RMSE_mis_P_PI, RMSE_cor_N_PI, RMSE_mis_N_PI, RMSE_cor_IF, RMSE_mis_IF)
  return(toreturn)
}




fun.simulate(200)
# testing = fun.simulate(1000)
# barplot(testing)
# 
# thelength = length(fun.simulate(1000))
# 
# # Repeat simulation for different sample sizes
samplesizes = c(rep(200, 30), rep(1000, 30), rep(5000, 30))
samplesizes = c(rep(200, 10))

arr <- array(dim=c(length(samplesizes),thelength+1))
colnames(arr) <- c("sample_sizes", 
                   "RMSE_cor_P_PI", "RMSE_mis_P_PI", "RMSE_cor_N_PI", 
                   "RMSE_mis_N_PI", "RMSE_cor_IF", "RMSE_mis_IF")
arr[1:length(samplesizes),1] = samplesizes
arr


for(s in 1:length(samplesizes)){
  fun.results = fun.simulate(samplesizes[s])
  arr[s, 2] = fun.results[1]
  arr[s, 3] = fun.results[2]
  arr[s, 4] = fun.results[3]
  arr[s, 5] = fun.results[4]
  arr[s, 6] = fun.results[5]
  arr[s, 7] = fun.results[6]
}
arr
df.sim = as.data.frame(arr)
arr1=arr



s=1
for(s in 1:length(samplesizes)){
    fun.results = fun.simulate(samplesizes[s])
    fun.results
    errorcount = which(is.na(fun.results[5]))
    
    dff[which(dff$A==0),], 
    
    arr[s, 2] = fun.results[1]
    arr[s, 3] = fun.results[2]
    arr[s, 4] = fun.results[3]
    arr[s, 5] = fun.results[4]
    
    while(!is.na(fun.results[5]) && !is.na(fun.results[6])){
    
    arr[s, 6] = fun.results[5]
    arr[s, 7] = fun.results[6]
    
    i=i+1
    }
  }

arr
arr
df.sim = as.data.frame(arr)
df.sim1 = as.data.frame(arr)
arr1=arr

colnames(df.sim)[2] = "Parametric_correctly_specified"
colnames(df.sim)[3] = "Parametric_misspecified"
colnames(df.sim)[4] = "Nonparametric_untransformed_Xs"
colnames(df.sim)[5] = "Nonparametric_transformed_Xs"
colnames(df.sim)[6] = "Influence_function_untransformed_Xs"
colnames(df.sim)[7] = "Influence_function_transformed_Xs"
# Backup table: 
# df.sim = read.table("/Users/mariacuellar/Desktop/table1.csv", sep=",", header=TRUE)[,2:5]
# head(df.sim)

setwd(WD_simulation)
save(df.sim, file="df.sim.fin-2017-05-15.rda")


# Alert when it's done running
beep(1) 


# Plotting the RMSE to compare estimators
mdf = melt(df.sim, id="sample_sizes")
mdf$Estimator = mdf$variable

setwd(WD_figs)
pdf(file="20170515_RMSE_comparison.pdf", width=10, height=5)
ggplot(na.omit(mdf), aes(x = sample_sizes, y=value, color=Estimator)) + geom_point() + 
  geom_smooth(method = "loess", size = 1.5) +
  labs(title = "Comparison of four estimators, 100 iterations per sample size", y="Root mean squared error", x="Sample size")
dev.off()



# Barplot with error bars
# Reshape data frame from simulatin, df.sim
aa = melt(data = df.sim, id.vars = "sample_sizes", 
          variable.name = "Algorithm", value.name = "RMSE") # Not changing the names right...
dat = rename(aa, c(sample_sizes = "Sample_Size", variable = "Algorithm", value = "RMSE"))
dat$Sample_Size = as.factor(dat$Sample_Size)

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

# Barplot with ggplot
bp = ggplot(dfsum, aes(x=Sample_Size, y=RMSE, fill=Algorithm)) + 
  geom_bar(stat="identity", color="black", position=position_dodge()) +
  geom_errorbar(aes(ymin=RMSE-sd, ymax=RMSE+sd), width=.2, position=position_dodge(.9)) +
  ggtitle("Comparison of estimators, 100 iterations per sample size") + xlab("Sample size") + ylab("Rootmean squared error")


setwd(WD_figs)
pdf("20170515_Barplot.pdf", width=8, height=3.3)
bp + scale_fill_brewer(palette="Paired")
dev.off()

# Nice colors
bp + scale_fill_manual(values=wes_palette(n=3, name="GrandBudapest"))
bp + scale_fill_brewer(palette="Paired")
.






# FITTING OTHER NONPARAMETRIC PLUG IN ESTIMATORS



#   
#   # Correctly specified model: NONparametric plug-in (P_PI) estimator:
#   # Kernel----
#   # Fitting model, a bit faster than without tol, ftol, but still slow
#   # got bandwidths from running it once on a sample of the data
#   cor_N_PI_mu0k = npreg(Y~X1+X2+X3+X4, data = dff[which(dff$A==0),], tol=0.1, ftol=0.1, 
#                         bws=c(25156071.0125331692, 0.6231756800, 0.9547994939, 0.2197950283))
#   
#   cor_N_PI_mu1k = npreg(Y~X1+X2+X3+X4, data = dff[which(dff$A==1),], tol=0.1, ftol=0.1, 
#                         bws=c(1063281.762, 2113100.136, 2943341.821, 1139805.992))
#   
#   # Getting fitted values
#   cor_N_PI_mu0_hatk = predict(cor_N_PI_mu0k, newdata = dff)#$mean # with faster kernel
#   cor_N_PI_mu1_hatk = predict(cor_N_PI_mu1k, newdata = dff)#$mean
#   
#   # Gamma hat
#   cor_N_PI_gammahatk = (cor_N_PI_mu1_hatk - cor_N_PI_mu0_hatk)/cor_N_PI_mu1_hatk
#   
#   # RMSE 
#   RMSE_cor_N_PIk = sqrt(  mean( (cor_N_PI_gammahatk - gamma)^2 )  ) # 0.2458764
#   RMSE_cor_N_PI = RMSE_cor_N_PIk



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



#   # Ranger---- 
#   # Fitting model
#   mis_N_PI_mu0r = ranger(Y~X1star+X2star+X3star+X4star, data = dff[which(dff$A==0),])
#   mis_N_PI_mu1r = ranger(Y~X1star+X2star+X3star+X4star, data = dff[which(dff$A==1),])
#   
#   # Getting fitted values (after inverse link function)  
#   mis_N_PI_mu0_hatr = expit(predict(mis_N_PI_mu0r, data=dff)[1][[1]])
#   mis_N_PI_mu1_hatr = expit(predict(mis_N_PI_mu1r, data=dff)[1][[1]])
#   
#   # Gamma hat
#   mis_N_PI_gammahatr = (mis_N_PI_mu1_hatr - mis_N_PI_mu0_hatr)/mis_N_PI_mu1_hatr
#   
#   # RMSE
#   RMSE_mis_N_PIr = sqrt(  mean( (mis_N_PI_gammahatr - gamma)^2 )  )
#   




# SuperLearner----
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
# mis_N_PI_mu0_hatSL = predict(mis_N_PI_mu0SL, newdata=as.data.frame(cbind(dff$X1star,dff$X2star,dff$X3star,dff$X4star)))$pred
# mis_N_PI_mu1_hatSL = predict(mis_N_PI_mu1SL, newdata=as.data.frame(cbind(dff$X1star,dff$X2star,dff$X3star,dff$X4star)))$pred
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
# Support vector machine----
# Fitting model:
cor_N_PI_mu0svm = svm(Y~X1+X2+X3+X4, data = dff[which(dff$A==0),])
cor_N_PI_mu1svm = svm(Y~X1+X2+X3+X4, data = dff[which(dff$A==1),])

# Getting predicted values
cor_N_PI_mu0_hatsvm = expit(predict(cor_N_PI_mu0svm, newdata = dff))
cor_N_PI_mu1_hatsvm = expit(predict(cor_N_PI_mu1svm, newdata = dff))

# Gamma hat
cor_N_PI_gammahatsvm = (cor_N_PI_mu1_hatsvm - cor_N_PI_mu0_hatsvm)/cor_N_PI_mu1_hatsvm

# RMSE
RMSE_cor_N_PIsvm = sqrt(  mean( (cor_N_PI_gammahatsvm - gamma)^2 )  ) # 0.3977806
RMSE_cor_N_PI = RMSE_cor_N_PIsvm


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
