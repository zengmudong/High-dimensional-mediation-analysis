library(glmnet)
library(stringr)
library(MASS)
library(scalreg)
library(Matrix)
library(freebird)
library(tidyr)
library(ggplot2)
library(dplyr)
library(latex2exp)
library(ggpubr)
source("utils_mediation.R")
# Simulation Studies:
## Mediation model with confounding variables
## Compare our proposed new method, oracle and ZWZ's method estimate \hat{\sigma}_1
## Paper section 4.2 Figure 3

## Mediation model with confounding variables
## Paper section 4.2
# simulation cases: 
#    - cases =  1:21, holding c2 = 0.5, c1= -1.0, -0.9, ..., -0.1, 0, 0.1, ..., 0.9, 1.0 
#    - cases = 22:42, holding c1 = 0.5, c2= -1.0, -0.9, ..., -0.1, 0, 0.1, ..., 0.9, 1.0 
cases= 1:42
# number of replications in the simulation
nsim = 500 

gen_error<-function(N,p,rho){
  X = matrix(NA,N,p)
  X[,1] = rnorm(N)
  for(ii in 2:p){
    X[,ii] = rho*X[,(ii-1)] + sqrt(1-rho^2)*rnorm(N)
  }
  return(X)
} 

simulations <- function(kk){
  #if (kk%%100 ==0){print(kk)} 
  alpha1 = c2 
  Gamma = coefs$Gamma
  Gamma[1,] = c1 * Gamma[1,]
  
  M = W%*%Gamma +  gen_error(n,p,rho = 0.5)
  M= scale(M)
  Y = M %*% alpha0 + X%*% alpha1 + 
    S%*% alpha2 + rnorm(n,0, sigma1) # change here 
  
  gamma_hat = solve(t(X)%*%X)%*%t(X)%*%Y
  Sigma_XX = t(X)%*%X/n
  invXX = solve(Sigma_XX) 
  Our_start = Sys.time()
  
  hbic= c() #matrix(NA, ncol = ngrid,nrow = length(rho2_grid))
  alpha0_rcd = matrix(NA,nrow = ngrid,ncol=p)
  sigma1_rcd=c()
  for( ii in 1: ngrid){
    result = HBIC_calc(lamb_grid[ii], X,Y,M, S =S, n_imp = 8 )
    hbic[ii] = result$BIC
    alpha0_rcd[ii,]=result$alpha0
    sigma1_rcd[ii] = result$sigma1_hat
    #print(which(result$alpha0[1:1000]!=0))
  }
  
  id = which(hbic==min(hbic))
  #plot(lamb_grid, hbic)
  id = tail(id,1)
  
  alpha0_hat = alpha0_rcd[id,]
  
  ## Use LLA to select model then refit using OLS
  A = which(alpha0_hat!=0)
  
  M_A = M[,A]
  Refit = mediationInference(X, Y, M_A, S) 
  # Oracle 
  Ocl= summary(lm(Y~0+cbind(X,M[,which(coefs$alpha0!=0)], S)))
  # Zhou et al. (2020)'s method
  result_scalreg = scalreg(cbind(M,X,S),Y)
  ZWZ_alpha_hat = result_scalreg$co
  ZWZ_sigma1_hat = result_scalreg$hsigma
  
  return(list(c1 = c1, c2= c2, R_sigma1 = Refit$sigma1_hat, Ocl_sigma1 = Ocl$sigma, ZWZ_sigma1 = ZWZ_sigma1_hat))
}
if(!file.exists("results/simulation_allS.Rdata")){
  print("Error: no exist simulation_allS.Rdata in results folder. Please run simulationPreprocessing.R before running this file.")
}
load("results/simulation_allS.Rdata")
alpha0 = coefs$alpha0
# set the true value for alpha0
alpha0[which(alpha0!=0)] = c(1.2,1.1,1.0,
                               -0.9,-0.8,-0.7,0.6,0.5,0.4,0.3,0.2)
A_true = which(alpha0!=0)
alpha2 = coefs$alpha2
sigma1 = coefs$sigma1
Gamma = coefs$Gamma

X = data.matrix(dt$X)
S= dt$S
W = cbind(X, S)
M = dt$M
n = nrow(W)
p = ncol(dt$M)
q = ncol(X)
s = ncol(dt$S)

beta0 = t(alpha0)%*%Gamma["trauma",]
# Setup HBIC for different 
ngrid = 20 # ngrid for lambda
lamb_grid = seq(0.13,0.25,length.out = ngrid)
lamb_grid0 = seq(0.1,0.3,length.out = ngrid)

c_gr = matrix(c(seq(-1, 1, 0.1), rep(0.5, 42), seq(-1,1,0.1)),ncol=2)
nc = nrow(c_gr)


for( jj in cases){
  c1 = c_gr[jj,1]
  c2 = c_gr[jj,2]
  cat("Running case ", jj, ', c1=', c1, ', c2 =',c2, '\n')
  
  tab = lapply(1:nsim, simulations)
  df = data.frame(matrix(unlist(tab),ncol = length(tab[[1]]),byrow = T))
  names(df) = names(tab[[1]])
  if (jj>1){
    df_old = rbind(df_old, df)
  }else{
    df_old = df
  }
  write.csv(df_old, "results/sigma1_compare.csv", row.names = FALSE,col.names = TRUE)
}

dt = read.csv("results/sigma1_compare.csv")
colnames(dt)<- c("c1", "c2", "New.method","Oracle", "ZWZ")
data_long <- pivot_longer(dt[1:(21*nsim),], cols =New.method:ZWZ, names_to="Methods", values_to = "sigma1_hat")
data_long$Methods = factor(data_long$Methods,levels = c("New.method","Oracle", "ZWZ"))
data_long$c1 = factor(data_long$c1)
Sn<-ggplot(data = data_long, aes(x = c1, y =sigma1_hat,fill=Methods)) + 
  geom_boxplot(outlier.alpha = 0, lwd=0.3)+
  theme_bw() +theme(panel.grid =element_blank()) +
  #geom_point(shape=16,position=position_jitterdodge(),size =0.2,alpha=0.5)+
  labs(y=TeX("$\\hat{\\sigma}_1$"), x= TeX("$c_1$"))+
  scale_x_discrete(breaks = as.character(seq(-1,1,0.2)),drop = F)+  
  scale_fill_manual(values=c("chartreuse3", "coral2", "cornflowerblue"))

Sn
data_long <- pivot_longer(dt[(21*nsim+1):nrow(dt),], cols =New.method:ZWZ, names_to="Methods", values_to = "sigma1_hat")
data_long$Methods = factor(data_long$Methods,levels = c("New.method","Oracle", "ZWZ"))
data_long$c2 = factor(data_long$c2)

Tn<-ggplot(data = data_long, aes(x = c2, y =sigma1_hat,fill=Methods)) + 
  geom_boxplot(outlier.alpha = 0, lwd=0.3)+
  theme_bw() +theme(panel.grid =element_blank()) +
  labs(y=TeX("$\\hat{\\sigma}_1$"), x= TeX("$c_2$"))+
  scale_x_discrete(breaks = as.character(seq(-1,1,0.2)),drop = F) +
  scale_fill_manual(values=c("chartreuse3", "coral2", "cornflowerblue"))
Tn
ggarrange(Sn, Tn, ncol = 2, 
          common.legend = TRUE, legend = "bottom")
ggsave("sigma1.pdf", dpi = 600, 
       width = 8, height =3,units = 'in')

