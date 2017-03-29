# Thesis simulations
# March 24, 2017
# Making my own simulation for the parameter of interest: PC = gamma = 1-1/RR = 1-mu0/mu1.

# define working directory
WD_figs = "/home/mcuellar"
WD_thesis = "/home/mcuellar/Thesis"

# install.packages("np")
# install.packages("beepr")

library(clusterPower) # lets me use expit function
library(dplyr) # lets me use sample_n function
library(hydroGOF) # lets me use rmse function
library(np) # lets me use nonparametric models, like kernel (https://cran.r-project.org/web/packages/np/np.pdf)
library(randomForest) # lets me do random forests

# Function to simulate estimation of PC (gamma) as plugin, nonparametric, parametric IF, and nonparametric IF

# samplesize = 1000 # test
fun.simulate = function(samplesize){
  # seed for random number generator
  #set.seed(400)
  
  # true parameters
  #samplesize = 1000
  index = 1:samplesize
  beta = 0.5
  Y1 = rbinom(n = samplesize, size = 1, prob = beta) #rbinom(n = samplesize, size = 1, prob = expit(210 + 27.4*X1 + 13.7*X2 + 13.7*X3 + 13.7*X4)) 
  
  X1 = rnorm(n = samplesize, mean = 0, sd = 1) # correct model
  X2 = rnorm(n = samplesize, mean = 0, sd = 1)
  X3 = rnorm(n = samplesize, mean = 0, sd = 1)
  X4 = rnorm(n = samplesize, mean = 0, sd = 1)
  
  X1star = exp(X1/2) # misspecified model
  X2star = X2/(1 + exp(X1)) + 10
  X3star = (X1*X3/25 + 0.6)^3
  X4star = (X2 + X4 + 20)^2
  
  A = rbinom(n = samplesize, size = 1, prob = expit(-X1+0.5*X2-0.25*X3-0.1*X4))
  gamma = expit(-X1+0.5*X2-0.25*X3-0.1*X4)
  
  # generate data frame
  df = as.data.frame(cbind(index, X1, X2, X3, X4, X1star, X2star, X3star, X4star, gamma, A, Y1))
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
  
  # ordering data so it's as it was at the beginning
  dff = df_wy0[order(df_wy0$index),] 
  head(dff)
  
  
  
  # now getting into the models ---
  
  # Correctly specified model: parametric plug-in (P_PI) estimator:
  cor_P_PI_mu0 = glm(Y~X1+X2+X3+X4, data = dff[which(dff$A==0),], family = "binomial") #fit only on items with A=0
  cor_P_PI_mu1 = glm(Y~X1+X2+X3+X4, data = dff[which(dff$A==1),], family = "binomial") #fit only on items with A=1
  
  # getting fitted values (after inverse link function)
  cor_P_PI_mu0_hat = expit(predict.glm(cor_P_PI_mu0, newdata = dff))
  cor_P_PI_mu1_hat = expit(predict.glm(cor_P_PI_mu1, newdata = dff))
  
  # calculating gamma hat
  cor_P_PI_gammahat = (cor_P_PI_mu1_hat - cor_P_PI_mu0_hat)/cor_P_PI_mu1_hat
  
  # RMSE 
  RMSE_cor_P_PI = sqrt(  mean( (cor_P_PI_gammahat - gamma)^2 )  )


  
  # Misspecified model: parametric plug-in (P_PI) estimator:
  mis_P_PI_mu0 = glm(Y~X1star+X2star+X3star+X4star, data = dff[which(dff$A==0),], family = "binomial") #fit only on items with A=0
  mis_P_PI_mu1 = glm(Y~X1star+X2star+X3star+X4star, data = dff[which(dff$A==1),], family = "binomial") #fit only on items with A=1
  
  # getting fitted values (after inverse link function)
  mis_P_PI_mu0_hat = expit(predict.glm(mis_P_PI_mu0, newdata = dff))
  mis_P_PI_mu1_hat = expit(predict.glm(mis_P_PI_mu1, newdata = dff))

  # calculating gamma hat
  mis_P_PI_gammahat = (mis_P_PI_mu1_hat - mis_P_PI_mu0_hat)/mis_P_PI_mu1_hat
  
  # RMSE 
  RMSE_mis_P_PI = sqrt(  mean( (mis_P_PI_gammahat - gamma)^2 )  )
  

  
  
  # Correctly specified model: NONparametric plug-in (P_PI) estimator:
  cor_N_PI_mu0 = randomForest(factor(Y)~X1+X2+X3+X4, data = dff[which(dff$A==0),])
  cor_N_PI_mu1 = randomForest(factor(Y)~X1+X2+X3+X4, data = dff[which(dff$A==1),])
  
  # getting fitted values (after inverse link function)
  cor_N_PI_mu0_hat = as.numeric(predict(cor_N_PI_mu0, newdata=dff))
  cor_N_PI_mu1_hat = as.numeric(predict(cor_N_PI_mu1, newdata=dff))
  
  # calculating gamma hat
  cor_N_PI_gammahat = (cor_N_PI_mu1_hat - cor_N_PI_mu0_hat)/cor_N_PI_mu1_hat
  
  # RMSE 
  RMSE_cor_N_PI = sqrt(  mean( (cor_N_PI_gammahat - gamma)^2 )  )
  
  
  
  # Misspecified model: NONparametric plug-in (P_PI) estimator:
  mis_N_PI_mu0 = randomForest(factor(Y)~X1star+X2star+X3star+X4star, data = dff[which(dff$A==0),])
  mis_N_PI_mu1 = randomForest(factor(Y)~X1star+X2star+X3star+X4star, data = dff[which(dff$A==1),])
  
  # getting fitted values (after inverse link function)
  mis_N_PI_mu0_hat = as.numeric(predict(cor_N_PI_mu0, newdata=dff))
  mis_N_PI_mu1_hat = as.numeric(predict(cor_N_PI_mu1, newdata=dff))
  
  # calculating gamma hat
  mis_N_PI_gammahat = (mis_N_PI_mu1_hat - mis_N_PI_mu0_hat)/mis_N_PI_mu1_hat
  
  # RMSE 
  RMSE_mis_N_PI = sqrt(  mean( (mis_N_PI_gammahat - gamma)^2 )  )
  

  
  # Function will return this
  return(c(RMSE_cor_P_PI, RMSE_mis_P_PI, RMSE_cor_N_PI, RMSE_mis_N_PI))
}

fun.simulate(1000)

# repeat simulation for different n's
#samplesizes = c(500, 1000, 1500, 2000, 2500, 3000, 3500, 10000)
samplesizes = c(500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 
                1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000,
                1500, 1500, 1500, 1500, 1500, 1500, 1500, 1500, 1500, 1500)

samplesizestot = seq(from = 500,to = 1500,by = 500)

arr <- array(dim=c(length(samplesizes),5))
colnames(arr) <- c("sample_sizes", "cor_P_PI_RMSE", "mis_P_PI_RMSE", "cor_N_PI_RMSE", "mis_N_PI_RMSE")
arr[1:length(samplesizes),1] = samplesizes
arr

for(s in 1:length(samplesizes)){
  arr[s, 2] = fun.simulate(samplesizes[s])[1]
  arr[s, 3] = fun.simulate(samplesizes[s])[2]
  arr[s, 4] = fun.simulate(samplesizes[s])[3]
  arr[s, 5] = fun.simulate(samplesizes[s])[4]
}

df.sim = as.data.frame(arr)
df.sim

setwd(WD_thesis)
save(df.sim, file="df.sim.fin-2017-03-29.rda")

# Making histograms
dat_n500 = df.sim[which(df.sim$sample_sizes==500),]

# How to make these??
mean(dat_n500$cor_P_PI_RMSE)
mean(dat_n500$mis_P_PI_RMSE)

mean(dat_n500$cor_N_PI_RMSE)
mean(dat_n500$mis_N_PI_RMSE)


# Scatterplot
setwd(WD_thesis)
pdf("MSE_estimating_PC_samesamplesize1.pdf", width=8, height=6)
par(mfrow=c(1,1))
plot(df.sim$sample_sizes, df.sim$P_PI_simulations, col="blue", lwd=2, ylim=c(0,0.06), 
     main="Root Mean Squared Error (RMSE) of PC Estimation", ylab="Mean Squared Error (MSE)", xlab="Sample size")
points(df.sim$sample_sizes, df.sim$N_PI_simulations, col="red", lwd=2)
panel.first = grid(lwd=1, col="gray30")
legend(x = 6000, y = 0.01, cbind("Parametric PI estimates", "Nonparametric PI estimates"), cex=0.8, lwd=2, col=c("blue", "red"),bg = "white")
dev.off()

df.sim.fin
df.sim.fin[17,3]=df.sim.fin[16,3]
df.sim.fin[34,3]=df.sim.fin[33,3]
df.sim.fin[39,3]=df.sim.fin[38,3]

# Boxplots

ddf = df.sim.fin[which(df.sim.fin$sample_sizes==samplesizestot[1]),]

samplesizestot = c(1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000, 10000)
arr_P <- array(dim=c(10,10))
colnames(arr_P) <- samplesizestot
arr_P

for(s in 1:length(samplesizestot)){
  arr_P[1:length(samplesizestot), s] = df.sim.fin[which(df.sim.fin$sample_sizes==samplesizestot[s]),]$P_PI_simulations
}

arr_N <- array(dim=c(10,10))
colnames(arr_N) <- samplesizestot
arr_N

for(s in 1:length(samplesizestot)){
  arr_N[1:length(samplesizestot), s] = df.sim.fin[which(df.sim.fin$sample_sizes==samplesizestot[s]),]$N_PI_simulations
}

pdf(file="Boxplots.pdf", width=10, height=5)
par(mfrow=c(1,2))
boxplot(arr_P, col = "yellow", ylim=c(0, 0.05), main="Parametric PI")
boxplot(arr_N, col = "yellow", ylim=c(0, 0.05), main="Nonparametric PI")
dev.off()


