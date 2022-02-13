library(glmnet)
library(stringr)
library(MASS)
library(Matrix)
library(globaltest)
library(scalreg)
library(freebird)
source("utils_mediation.R")
# Simulation Studies:
## Mediation model with confounding variables
## Paper section S.3
# Set grid points for c2
c2_grid = c(0.5, 0.25, 0.125, 0.0625)
# simulation cases: 
#    - cases =  1:21, holding c2, c1= -1.0, -0.9, ..., -0.1, 0, 0.1, ..., 0.9, 1.0 (Test for indirect effect Sn)
cases= 1:21
# number of replications in the simulation
nsim = 500


gen_error<-function(N,p,rho){
  ## Generate AR(1) Error matrix, (i,j)th element is \rho^{|i-j|}
  X = matrix(NA,N,p)
  X[,1] = rnorm(N)
  for(ii in 2:p){
    X[,ii] = rho*X[,(ii-1)] + sqrt(1-rho^2)*rnorm(N)
  }
  return(X)
} 


mediationInference3<-function(X, Y, M, S = NULL, M_tld = NULL, 
                    alpha0_hat = NULL, alpha0_tld = NULL, alpha1_hat = NULL,alpha2_hat=NULL){
  p = ncol(M)
  q = ncol(X)
  n = nrow(X)
  if(is.null(S)){
    Z = cbind(M,X)
    s = 0
    if(is.null(M_tld)){MS = M }else{ MS = M_tld}
    if(is.null(alpha0_tld)){alpha0_tld = solve(t(MS)%*%MS)%*%t(MS)%*%Y}
    RSS02 = t(Y - MS%*% alpha0_tld) %*% (Y - MS%*%alpha0_tld)
  }
  else{
    Z = cbind(M,X,S)
    s = ncol(S)
    if(is.null(M_tld)){MS = cbind(M,S) }else{ MS = cbind(M_tld,S)}
    if(is.null(alpha0_tld)){alpha0_tld = solve(t(MS)%*%MS)%*%t(MS)%*%Y}
    RSS02 = t(Y - MS%*% alpha0_tld) %*% (Y - MS%*%alpha0_tld)
  }
  if(is.null(alpha0_hat)){
    alpha_rf = solve(t(Z)%*%Z)%*%t(Z)%*%Y
    alpha0_hat = alpha_rf[1:p]
    alpha1_hat = alpha_rf[(p+1):(p+q)]
    if(s >0){alpha2_hat = alpha_rf[(p+q+1):(p+q+s)]}
  }
  else{
    alpha_rf=c(alpha0_hat,alpha1_hat,alpha2_hat)
  }
  
  res = Y - Z%*%alpha_rf
  RSS12 = as.numeric(t(res) %*% (res)) # Test direct effect
  
  # RSS01 = t(Y - X%*%alpha1_hat) %*% (Y - X%*%alpha1_hat) # Test indirect effect
  # RSS11 = t(Y-X%*%gamma_hat) %*% (Y-X%*%gamma_hat)
  
  df = p+q + s
  sigma1_hat = RSS12/(n - df)
  Sigma_MM = t(M)%*%M /n
  if(s == 0){
    gamma_hat = solve(t(X)%*%X)%*%t(X)%*%Y 
    sigmaT_hat = t(Y-X%*%gamma_hat) %*% (Y-X%*%gamma_hat)/(n-q)
    beta_hat = gamma_hat -alpha1_hat
    sigma2_hat = pmax(0,(sigmaT_hat - sigma1_hat))
    
    #tmp1 = cbind(t(X)%*%X,t(X)%*%M)
    #tmp2 = cbind(t(M)%*%X,t(M)%*%M )
    #Sigma_hat = rbind(tmp1,tmp2)/n
    invXX = solve(t(X)%*%X/n)
    Sigma_MX =t(M)%*%X/n
    
    B = invXX %*%t(Sigma_MX) %*%solve(Sigma_MM -  Sigma_MX%*%invXX %*% t(Sigma_MX)) %*%Sigma_MX %*%invXX 
    var_alpha1_hat = sigma1_hat*(invXX + B)
    cov_beta_hat = sigma2_hat * invXX + sigma1_hat * B
  }else{
    V = cbind(X,S)
    gamma_hat = solve(t(V)%*%V)%*%t(V)%*%Y 
    sigmaT_hat = t(Y-V%*%gamma_hat) %*% (Y-V%*%gamma_hat)/(n-q-s)
    beta_hat = gamma_hat[1:q] -alpha1_hat
    sigma2_hat = pmax(0,(sigmaT_hat - sigma1_hat))
    Sigma_VV = t(V)%*%V/n
    invVV = solve(Sigma_VV)
    Sigma_MV =t(M)%*%V/n
    Sigma_VM =t(Sigma_MV)
    B = invVV %*% Sigma_VM %*%solve(Sigma_MM - Sigma_MV %*% invVV %*% Sigma_VM)%*%Sigma_MV%*%invVV
    var_alpha1_hat = sigma1_hat*(invVV + B)[1:q,1:q]
    cov_beta_hat = sigma2_hat * invVV + sigma1_hat * B
    cov_beta_hat = cov_beta_hat[1:q, 1:q]
  }
  # Test for beta
  # Wald's test
  Sn = n*t(beta_hat) %*% solve(cov_beta_hat) %*% beta_hat
  # Test for alpha1
  # LRT
  Tn1 = (n-df) * (RSS02-RSS12)/RSS12
  #Tn2 = n*log(RSS02/RSS12)
  return(list(Sn = Sn, Tn = Tn1, 
              beta_hat = beta_hat, alpha0_hat = alpha0_hat,
              alpha1_hat = alpha1_hat, alpha2_hat = alpha2_hat, B = B, 
              var_beta = cov_beta_hat, var_alpha1_hat = var_alpha1_hat))
}


Vera_pcd<- function(X,Y,M, S2, alpha2, Gamma){
  colnames(M) <- paste0("M", 1:ncol(M))
  n = dim(M)[1]
  # Null model contain intercept and X
  # method2: oracle
  res_Y = Y - S2 %*% alpha2
  res_M = M - S2%*% Gamma[3:nrow(Gamma),]
  MY_test <- gt(response= res_Y,  
                alternative =res_M , 
                null = ~X, model = 'linear')
  # Null model contain intercept
  XM_test <- gt(X,res_M, model ='linear') # M = X + S
  H.pvalue = pmax(MY_test@result[1], XM_test@result[1])

  return(H.pvalue)
}

simulations <- function(kk){
  
  set.seed(as.numeric(as.POSIXct(Sys.time())))
  # Generate data
  alpha1 = c2 
  #Gamma = c1 * coefs$Gamma[1,]
  #M = X%*%Gamma +  gen_error(n,p,rho = 0.5)
  #Y = M %*% alpha0 + X%*% alpha1 + rnorm(n,0, sigma1) 
  
  Gamma = coefs$Gamma
  Gamma[1,] = c1 * Gamma[1,]
  M = W%*%Gamma +  gen_error(n,p,rho = 0.5)
  Y = M %*% alpha0 + X%*% alpha1 +    S%*% alpha2 + rnorm(n,0, sigma1)
  
  gamma_hat = solve(t(X)%*%X)%*%t(X)%*%Y
  Sigma_XX = t(X)%*%X/n
  invXX = solve(Sigma_XX) 
  Our_start = Sys.time()
  # HBIC to select best tuning parameters
  hbic= c() 
  results =lapply(lamb_grid, HBIC_calc, xx=X,yy=Y,mm=M,S=S,n_imp = 8) ##
  for( ii in 1: ngrid){
    hbic[ii] = results[[ii]]$BIC
  }
  
  id = which(hbic==min(hbic))
  #plot(lamb_grid, hbic)
  id = tail(id,1)
  result = results[[id]]
  alpha0_hat = result$alpha0
  alpha1_hat = result$alpha1
  intcpt = matrix(rep(1, n), ncol = 1)
  hbic0= c() 
  results =lapply(lamb_grid0, HBIC_calc, xx=S ,yy=Y,mm=M,n_imp = 8) ##
  
  for( ii in 1: ngrid){
    hbic0[ii] = results[[ii]]$BIC
  }
  
  id = which(hbic0==min(hbic0))  
  
  #plot(lamb_grid, hbic)
  id = tail(id,1)
  alpha0_tld = results[[id]]$alpha0
  
  OurMethod_time = unclass(difftime(Sys.time(),Our_start, units = "secs"))[1]
  
  ## Use LLA to select model then refit using OLS
  A = which(alpha0_hat!=0)
  A_tld = which(alpha0_tld!=0)
  
  M_A = M[,A]
  Refit = mediationInference3(X, Y, M_A, S = S, M_tld =M[,which(alpha0_tld!=0)]) 
  # Oracle
  Ocl = mediationInference3(X, Y, M[,which(coefs$alpha0!=0)], S=S)
  
  # Zhou et al. (2020)'s method
  ZWZ_start = Sys.time()
  # Tuning parameters lam_list should change base on different settings c1 or c2 to achieve the best result
  ZWZ = hilma(Y,scale(M),cbind(X,S[,2:9]),mediation_setting='incomplete', lam_list = c(sqrt(log(p)/n)/3,sqrt(log(p)/n)/6))
  ZWZ_time = unclass(difftime(Sys.time(),ZWZ_start, units = "secs"))[1]
  
  
  Vera = Vera_pcd(X,Y,M, S[,2:9], alpha2[2:9], Gamma)
  return( list(Sn_rf=Refit$Sn, beta_hat_rf = Refit$beta_hat, var_beta_rf =Refit$var_beta, 
               Tn1_rf = Refit$Tn,
               alpha1_hat_rf = Refit$alpha1_hat, var_alpha1_rf = Refit$var_alpha1_hat, OurMethod_time = OurMethod_time,
               OclSn = Ocl$Sn, OclTn1 = Ocl$Tn,
               Oclbhat = Ocl$beta_hat, Oclvar_beta = Ocl$var_beta,
               Oclalpha1_hat = Ocl$alpha1_hat, Oclvar_alpha1_hat = Ocl$var_alpha1_hat,
               #ZWZbhat = ZWZ$beta_hat, ZWZ_var_beta = ZWZ$sigma_beta_hat, ZWZ_Sn = ZWZ$teststat_beta,
               #ZWZ_alpha1 = ZWZ$alpha1_hat, ZWZ_var_alpha1 = ZWZ$sigma_alpha1_hat, ZWZ_Tn = ZWZ$teststat_alpha1,
               #ZWZ_Sn_pvalue = ZWZ$pvalue_beta_hat, ZWZ_Tn_pvalue = ZWZ$pvalue_alpha1_hat, ZWZ_time = ZWZ_time,
               ZWZbhat = ZWZ$beta_hat[1], ZWZ_var_beta = ZWZ$sigma_beta_hat[1,1], ZWZ_Sn = ZWZ$teststat_beta[1],
               ZWZ_alpha1 = ZWZ$alpha1_hat[1], ZWZ_var_alpha1 = ZWZ$sigma_alpha1_hat[1,1], ZWZ_Tn = ZWZ$teststat_alpha1[1],
               ZWZ_Sn_pvalue = ZWZ$pvalue_beta_hat[1], ZWZ_Tn_pvalue = ZWZ$pvalue_alpha1_hat[1], ZWZ_time = ZWZ_time,
               Vera_H_pvalue = Vera
  ))
}

load("results/simulation_allS.Rdata")
alpha0 = coefs$alpha0
alpha0[which(alpha0!=0)] = c(1.0,0.9,0.8,
                             -0.9,-0.8,-0.7,0.6,0.5,0.4,0.3,0.2)
alpha2 = coefs$alpha2
sigma1 = coefs$sigma1
Gamma = coefs$Gamma

X = data.matrix(dt$X)
S= dt$S
W = cbind(X, S)
M = dt$M
n = nrow(M)
p = ncol(dt$M)
q = ncol(X)
s = ncol(S)

beta0 = t(alpha0)%*%Gamma["trauma",]
# Setup HBIC for different 
ngrid = 20 # ngrid for lambda
# The range of tuning parameters, search the best HBIC on an equal space grid
# The grid range may be adjusted for different settings to achieve the best result
lamb_grid = seq(0.2,0.39,length.out = ngrid) # Under full model
lamb_grid0 = seq(0.27,0.5,length.out = ngrid) # Under reduce model H0: \alpha_1=0

cn = c('c1','c2',
       'R beta Bias','R beta sample std','R beta est std','R beta MSE','R Sn Size/Power', 'R beta std SE',
       'R alpha1 Bias','R alpha1 sample std','R alpha1 est std','R alpha1 MSE','R Tn Size/Power', 'R alpha std SE',
       'O beta Bias','O beta sample std','O beta est std','O beta MSE','O Sn Size/Power', 'O beta std SE',
       'O alpha1 Bias','O alpha1 sample std','O alpha1 est std','O alpha1 MSE','O Tn Size/Power', 'O alpha std SE',
       'Z beta Bias','Z beta sample std','Z beta est std','Z beta MSE','Z Sn Size/Power', 'Z beta std SE',
       'Z alpha1 Bias','Z alpha1 sample std','Z alpha1 est std','Z alpha1 MSE','Z Tn Size/Power', 'Z alpha std SE',
       'Vera'
)

## There are three loops:

#### 1st loop: Whether the mediation model contain confounders or not

###### 2nd loop: c2 choose value from c2_grid: 0.5, 0.25, 0.125, 0.0625

####### 3rd loop: c1 choose value from seq(-1.0, 1.0, 0.1) hold the c2 not change

for (hasConfounder in c(TRUE, FALSE)){# 1st loop
  
  for(c2 in c2_grid){# 2nd loop
    if (hasConfounder){
      print("Running mediation model with confounding variables:")
      summary_fn = paste("results/simulation_confounder_globaltest_c2_",c2,".csv",sep="")
    }else{
      print("Running mediation model without confounding variables:")
      summary_fn = paste("results/simulation_globaltest_c2_",c2,".csv",sep="")
    }
    
    c1_range = seq(-1.0, 1.0, 0.1)
    c_gr = matrix(c(c1_range, rep(c2,length(c1_range))),ncol=2)
    nc = nrow(c_gr)
    
    for(jj in cases){# 3rd loop
      c1 = c_gr[jj,1]
      c2 = c_gr[jj,2]
      cat("Running case ", jj, ', c1=', c1, ', c2 =',c2, '\n')
      if (hasConfounder){
        individual_fn = paste("results/global_test_noConfounder_c1_",c1,"_c2_",c2,".csv",sep="")}
      else{
        individual_fn = paste("results/global_test_Confounder_c1_",c1,"_c2_",c2,".csv",sep="")
      }
      tab = lapply(1:nsim, simulations)
      df = data.frame(matrix(unlist(tab),ncol = length(tab[[1]]),byrow = T))
      names(df) = names(tab[[1]])
      
      if (file.exists(summary_fn)){
        result = read.csv(summary_fn,header = TRUE)
      }else{
        result = matrix(NA,ncol = length(cn),nrow = nc)
        colnames(result) = cn}
      
      result[jj,1] = c1
      result[jj,2] = c2
      
      result[jj,3] = mean(df$beta_hat_rf)-c1*beta0 # Refit bias of estiamted coeffcient beta
      result[jj,4] = sqrt(var(df$beta_hat_rf)) # Refit sample std
      result[jj,5] = sqrt(mean(df$var_beta_rf)/n) # Refit std estimator 
      result[jj,6] = var(df$beta_hat_rf) + (mean(df$beta_hat_rf)-c1*beta0)^2 # Refit beta MSE
      result[jj,7] = mean(df$Sn_rf>qchisq(0.95,1)) # Refit Sn size or power
      result[jj,8] = sqrt(var(sqrt(df$var_beta_rf/n))) # Refit std of estimated std for beta
      
      result[jj,9] = mean(df$alpha1_hat_rf)-c2 # Refit bias of estiamted coeffcient alpha1
      result[jj,10] = sqrt(var(df$alpha1_hat_rf)) # Refit sample std
      result[jj,11] = sqrt(mean(df$var_alpha1_rf)/n) # Refit std estimator 
      result[jj,12] = var(df$alpha1_hat_rf) + (mean(df$alpha1_hat_rf)-c2)^2 # Refit alpha1 MSE 
      result[jj,13] =mean(df$Tn1_rf>qchisq(0.95,1)) # Refit Tn1 size or power
      result[jj,14] = sqrt(var(sqrt(df$var_alpha1_rf/n))) # Refit std of estimated std for alpha1
      
      result[jj,15] = mean(df$Oclbhat)-c1*beta0 # Oracle bias of estiamted coeffcient beta 
      result[jj,16] = sqrt(var(df$Oclbhat)) # Oracle sample std 
      result[jj,17] = sqrt(mean(df$Oclvar_beta)/n) # Oracle std estimator 
      result[jj,18] = var(df$Oclbhat) + (mean(df$Oclbhat)-c1*beta0)^2 # Oracle beta MSE
      result[jj,19] = mean(df$OclSn>qchisq(0.95,1)) # Oracle Sn size or power
      result[jj,20] = sqrt(var(sqrt(df$Oclvar_beta/n))) # Oracle std of estimated std for beta
      
      result[jj,21] = mean(df$Oclalpha1_hat)-c2 # Oracle bias of estiamted coeffcient alpha1 
      result[jj,22] = sqrt(var(df$Oclalpha1_hat)) # Oracle alpha1 sample std
      result[jj,23] = sqrt(mean(df$Oclvar_alpha1_hat)/n) # Oracle alpha1 std estimator
      result[jj,24] = var(df$Oclalpha1_hat) + (mean(df$Oclalpha1_hat)-c2)^2 # Oracle alpha1 MSE
      result[jj,25] = mean(df$OclTn1>qchisq(0.95,1)) # Oracle Tn1 size or power
      result[jj,26] = sqrt(var(sqrt(df$Oclvar_alpha1_hat/n))) # Oracle std of estimated std for alpha1
      
      result[jj,27] = mean(df$ZWZbhat) - c1*beta0 # Zhou bias of estiamted coeffcient beta
      result[jj,28] = sqrt(var(df$ZWZbhat)) # Zhou sample std
      result[jj,29] = sqrt(mean(df$ZWZ_var_beta)/n) # Zhou beta est std
      result[jj,30] = (mean(df$ZWZbhat) - c1*beta0)^2 + var(df$ZWZbhat) # Zhou MSE of beta
      result[jj,31] = mean(df$ZWZ_Sn_pvalue<0.05) # Zhou Sn size/power
      result[jj,32] = sqrt(var(sqrt(df$ZWZ_var_beta/n))) # Zhou std of estimated std for beta
      
      result[jj,33] = mean(df$ZWZ_alpha1)-c2 # Zhou bias of estiamted coeffcient alpha1
      result[jj,34] = sqrt(var(df$ZWZ_alpha1)) # Zhou sample std
      result[jj,35] = sqrt(mean(df$ZWZ_var_alpha1)/n) # Zhou std for alpha1
      result[jj,36] = (mean(df$ZWZ_alpha1)-c2)^2 + var(df$ZWZ_alpha1)
      result[jj,37] = mean(df$ZWZ_Tn_pvalue<0.05) # Zhou Tn size / power
      result[jj,38] = sqrt(var(sqrt(df$ZWZ_var_alpha1/n))) # ZWZ std of estimated std for alpha1
      result[jj,39] =  mean(df$Vera_H_pvalue<0.05)
      write.csv(result,summary_fn,row.names = FALSE,col.names = TRUE)
      #write.csv(df, individual_fn, row.names = FALSE,col.names = TRUE)
    }
    
  }
}
  



