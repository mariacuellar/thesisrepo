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
library(np) # lets me use nonparametric models (https://cran.r-project.org/web/packages/np/np.pdf)
library(beepr)


# Function to simulate estimation of PC (gamma) as plugin, nonparametric, parametric IF, and nonparametric IF

# samplesize = 1000 # test
fun.simulate = function(samplesize){
  # seed for random number generator
  #set.seed(400)
  
  # true parameters
  samplesize = 1000
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
  df_gendat = df_wy0[order(df_wy0$index),] 
  head(df_gendat)
  
  
  
  # now getting into the models ---
  
  # Correctly specified model: parametric plug-in (P_PI) estimator:
  corr_P_PI = glm(Y~X1+X2+X3+X4, data = df_gendat, family = "binomial")
  
  ######### Crap, how do I get gamma hat?? ##########
  
  df_gendat$fittedY = fitted(corr_P_PI)
  
  df_gendat_Y11 = df_gendat[which(df_gendat$Y1==1),]
  
  
  
  gammahat = fitted(corr_P_PI)
  RMSE_P_PI = sqrt(  mean( (gammahat - gamma)^2 )  )
  
  
  # gettign fitted values (after inverse link function)
  fit_P_PI = fitted(corr_P_PI)

  
  # making them the same length
  cut_fit_P_PI_m0 = fit_P_PI_m0[1:(min(length(fit_P_PI_m0), length(fit_P_PI_m1)))]
  cut_fit_P_PI_m1 = fit_P_PI_m1[1:(min(length(fit_P_PI_m0), length(fit_P_PI_m1)))]
  
  # estimating gamma hat
  cutgammahat_P_PI = 1-cut_fit_P_PI_m0/cut_fit_P_PI_m1
  
  # MSE 
  cutgamma = gamma[1:(length(cutgammahat_P_PI))]
  MSE_P_PI = mean(  (cutgammahat_P_PI - cutgamma)^2  )
  
  
  
  # Nonparametric plug-in (N_PI) estimator
  N_PI_m0 = npregbw(Y~X, bandwidth.compute = TRUE, bwmethod = "cv.aic", data = subset(df_all, A == 0))
  N_PI_m1 = npregbw(Y~X, bandwidth.compute = TRUE, bwmethod = "cv.aic", data = subset(df_all, A == 1))
  
  # gettign fitted values (after inverse link function)
  fit_N_PI_m0 = predict(N_PI_m0)$mean
  fit_N_PI_m1 = predict(N_PI_m1)$mean
  
  # making them the same length
  cut_fit_N_PI_m0 = fit_N_PI_m0[1:(min(length(fit_N_PI_m0), length(fit_N_PI_m1)))]
  cut_fit_N_PI_m1 = fit_N_PI_m1[1:(min(length(fit_N_PI_m0), length(fit_N_PI_m1)))]
  
  # estimating gamma hat
  cutgammahat_N_PI = 1-cut_fit_N_PI_m0/cut_fit_N_PI_m1
  
  # MSE 
  cutgamma = gamma[1:(length(cutgammahat_N_PI))]
  MSE_N_PI = mean(  (cutgammahat_N_PI - cutgamma)^2  )
  
  
  return(c(MSE_P_PI, MSE_N_PI))
}


# repeat simulation for different n's
#samplesizes = c(500, 1000, 1500, 2000, 2500, 3000, 3500, 10000)
samplesizes = c(1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000,
                1500, 1500, 1500, 1500, 1500, 1500, 1500, 1500, 1500, 1500,
                2000, 2000, 2000, 2000, 2000, 2000, 2000, 2000, 2000, 2000,
                2500, 2500, 2500, 2500, 2500, 2500, 2500, 2500, 2500, 2500,
                3000, 3000, 3000, 3000, 3000, 3000, 3000, 3000, 3000, 3000,
                3500, 3500, 3500, 3500, 3500, 3500, 3500, 3500, 3500, 3500,
                4000, 4000, 4000, 4000, 4000, 4000, 4000, 4000, 4000, 4000,
                4500, 4500, 4500, 4500, 4500, 4500, 4500, 4500, 4500, 4500,
                5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000,
                10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000)

  # Making the list of samplesizes
  # ss = data.frame()
  # ss1 = data.frame()
  # theseq = seq(from = 1000,to = 5000,by = 500)
  # for(i in 1:length(theseq)){
  #   ss1 = rep((theseq)[i], times=10)
  #   ss = c(ss, ss1)
  # }
  # length(ss)
  # 
  # for(ii in 91:100){
  #   ss[[ii]] = 10000
  # }
  # sss=as.vector(unlist(ss))
  
  
samplesizestot = seq(from = 1000,to = 5000,by = 500)

arr <- array(dim=c(length(samplesizes),3))
colnames(arr) <- c("sample_sizes", "P_PI_simulations", "N_PI_simulations")
arr[1:length(samplesizes),1] = samplesizes
arr

for(s in 1:length(samplesizes)){
  arr[s, 2] = fun.simulate(samplesizes[s])[1]
  arr[s, 3] = fun.simulate(samplesizes[s])[2]
}

df.sim1 = as.data.frame(arr)
df.sim1

df.sim.fin = rbind(df.sim, df.sim1)
df.sim.fin

df.sim = df.sim.fin

setwd(WD_thesis)
save(df.sim.fin, file="df.sim.fin.rda")


# Scatterplot
setwd(WD_thesis)
pdf("MSE_estimating_PC_samesamplesize1.pdf", width=8, height=6)
par(mfrow=c(1,1))
plot(df.sim$sample_sizes, df.sim$P_PI_simulations, col="blue", lwd=2, ylim=c(0,0.06), 
     main=" Squared Error (MSE) of PC Estimation", ylab="Mean Squared Error (MSE)", xlab="Sample size")
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


