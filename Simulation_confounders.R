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
## Paper section 4.2
# simulation cases: 
#    - cases =  1:21, holding c2 = 0.5, c1= -1.0, -0.9, ..., -0.1, 0, 0.1, ..., 0.9, 1.0 (Test for indirect effect Sn)
#    - cases = 22:42, holding c1 = 0.5, c2= -1.0, -0.9, ..., -0.1, 0, 0.1, ..., 0.9, 1.0 (Test for direct effect Tn)
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

mediationInference2<-function(X, Y, M, S = NULL, M_tld = NULL){
  p = ncol(M)
  q = ncol(X)
  n = nrow(X)
  if(length(S) ==0){
    Z = cbind(M,X)
    s = 0
    alpha0_tld = solve(t(M)%*%M)%*%t(M)%*%Y
    RSS02 = t(Y - M%*% alpha0_tld) %*% (Y - M%*%alpha0_tld)}
  else{
    Z = cbind(M,X,S)
    s = ncol(S)
    if(is.null(M_tld)){MS = cbind(M,S) }
    else{ MS = cbind(M_tld,S)}
    alpha0_tld = solve(t(MS)%*%MS)%*%t(MS)%*%Y
    RSS02 = t(Y - MS%*% alpha0_tld) %*% (Y - MS%*%alpha0_tld)
  }
  alpha_rf = solve(t(Z)%*%Z)%*%t(Z)%*%Y
  alpha0_hat = alpha_rf[1:p]
  alpha1_hat = alpha_rf[(p+1):(p+q)]
  alpha2_hat = alpha_rf[(p+q+1):(p+q+s)]
  
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
  Tn1 = (n- q -s) * (RSS02-RSS12)/RSS12
  #Tn2 = n*log(RSS02/RSS12)
  return(list(Sn = Sn, Tn = Tn1, 
              beta_hat = beta_hat, alpha0_hat = alpha0_hat,
              alpha1_hat = alpha1_hat, alpha2_hat = alpha2_hat, B = B, 
              var_beta = cov_beta_hat, var_alpha1_hat = var_alpha1_hat))
}

simulations <- function(kk){
  if (kk%%100 ==0){print(kk)} 
  set.seed(1+kk)
  # Generate data
  alpha1 = c2 
  Gamma = coefs$Gamma
  Gamma[1,] = c1 * Gamma[1,]
  M = W%*%Gamma +  gen_error(n,p,rho = 0.5)
  Y = M %*% alpha0 + X%*% alpha1 + 
    S%*% alpha2 + rnorm(n,0, sigma1) # change here 
  
  gamma_hat = solve(t(X)%*%X)%*%t(X)%*%Y
  Sigma_XX = t(X)%*%X/n
  invXX = solve(Sigma_XX) 
  Our_start = Sys.time()
  
  hbic= c() 
  results =lapply(lamb_grid, HBIC_calc, xx=X,yy=Y,mm=M, S=S, n_imp = 8)
  for( ii in 1: ngrid){
    hbic[ii] = results[[ii]]$BIC
  }
  
  id = which(hbic==min(hbic))
  id = tail(id,1)
  result = results[[id]]
  alpha0_hat = result$alpha0
  alpha1_hat = result$alpha1
  
  hbic0= c() 
  results =lapply(lamb_grid0, HBIC_calc, xx=S ,yy=Y,mm=M,n_imp = 8)
  for( ii in 1: ngrid){
    hbic0[ii] = results[[ii]]$BIC
  }
  
  id = which(hbic0==min(hbic0))
  id = tail(id,1)
  alpha0_tld = results[[id]]$alpha0
  
  OurMethod_time = unclass(difftime(Sys.time(),Our_start, units = "secs"))[1]
  
  ## Use LLA to select model then refit using OLS
  A = which(alpha0_hat!=0)
  A_tld = which(alpha0_tld!=0)
  M_A = M[,A]
  # Our proposed new method
  Refit = mediationInference2(X, Y, M_A, S = S, M_tld =M[,which(alpha0_tld!=0)]) 
  # Oracle
  Ocl = mediationInference2(X, Y, M[,which(coefs$alpha0!=0)], S=S)

  # Zhou et al. (2020)'s method
  ZWZ_start = Sys.time()
  # Tuning parameters lam_list should change base on different settings c1 or c2 to achieve the best result
  ZWZ = hilma(Y,scale(M),cbind(X,S[,2:9]),mediation_setting='incomplete', lam_list = c(sqrt(log(p)/n)/3,sqrt(log(p)/n)/6))
  ZWZ_time = unclass(difftime(Sys.time(),ZWZ_start, units = "secs"))[1]
  
  return( list(Sn_rf=Refit$Sn, beta_hat_rf = Refit$beta_hat, var_beta_rf =Refit$var_beta, 
               Tn1_rf = Refit$Tn,
               alpha1_hat_rf = Refit$alpha1_hat, var_alpha1_rf = Refit$var_alpha1_hat, OurMethod_time = OurMethod_time,
               OclSn = Ocl$Sn, OclTn1 = Ocl$Tn,
               Oclbhat = Ocl$beta_hat, Oclvar_beta = Ocl$var_beta,
               Oclalpha1_hat = Ocl$alpha1_hat, Oclvar_alpha1_hat = Ocl$var_alpha1_hat,
               ZWZbhat = ZWZ$beta_hat[1], ZWZ_var_beta = ZWZ$sigma_beta_hat[1,1], ZWZ_Sn = ZWZ$teststat_beta[1],
               ZWZ_alpha1 = ZWZ$alpha1_hat[1], ZWZ_var_alpha1 = ZWZ$sigma_alpha1_hat[1,1], ZWZ_Tn = ZWZ$teststat_alpha1[1],
               ZWZ_Sn_pvalue = ZWZ$pvalue_beta_hat[1], ZWZ_Tn_pvalue = ZWZ$pvalue_alpha1_hat[1], ZWZ_time = ZWZ_time
  ))
}
if(!file.exists("results/simulation_allS.Rdata")){
  print("Error: no exist simulation_allS.Rdata in results folder. Please run simulationPreprocessing.R before running this file.")
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
n = nrow(W)
p = ncol(dt$M)
q = ncol(X)
s = ncol(dt$S)
beta0 = t(alpha0)%*%Gamma["trauma",]
# Setup HBIC for different 
ngrid = 20 # ngrid for lambda
# The range of tuning parameters, search the best HBIC on an equal space grid
# The grid range may be adjusted for different settings to achieve the best result
lamb_grid = seq(0.2,0.39,length.out = ngrid) # Under full model
lamb_grid0 = seq(0.27,0.46,length.out = ngrid) # Under reduce model H0: \alpha_1=0

c_gr = matrix(c(seq(-1, 1, 0.1), rep(0.5, 42), seq(-1,1,0.1)),ncol=2)
nc = nrow(c_gr)

cn = c('c1','c2',
       'R beta Bias','R beta sample std','R beta est std','R beta MSE','R Sn Size/Power', 'R beta std SE',
       'R alpha1 Bias','R alpha1 sample std','R alpha1 est std','R alpha1 MSE','R Tn Size/Power', 'R alpha std SE',
       'O beta Bias','O beta sample std','O beta est std','O beta MSE','O Sn Size/Power', 'O beta std SE',
       'O alpha1 Bias','O alpha1 sample std','O alpha1 est std','O alpha1 MSE','O Tn Size/Power', 'O alpha std SE',
       'Z beta Bias','Z beta sample std','Z beta est std','Z beta MSE','Z Sn Size/Power', 'Z beta std SE',
       'Z alpha1 Bias','Z alpha1 sample std','Z alpha1 est std','Z alpha1 MSE','Z Tn Size/Power', 'Z alpha std SE'
       
)
# Save final summary result file name
summary_fn = "results/simulation_confounders.csv"
for( jj in cases){
  c1 = c_gr[jj,1]
  c2 = c_gr[jj,2]
  cat("Running case ", jj, ', c1=', c1, ', c2 =',c2, '\n')
  
  individual_fn = paste("Gene_allS10908_JASA_c1_",c1,"_c2_",c2,".csv",sep="")
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
    
    write.csv(result,summary_fn,row.names = FALSE,col.names = TRUE)
  #}
  write.csv(df, individual_fn, row.names = FALSE,col.names = TRUE)
}

#############################################################################################################
## Plot Figure 2 in the paper section 4.2
summary_fn = "results/simulation_confounders.csv"

df = read.csv(summary_fn, header = T)
dt1 = df[1:21,c("c1","R.Sn.Size.Power", "O.Sn.Size.Power","Z.Sn.Size.Power")]
colnames(dt1) = c("c1","New.method", "Oracle", "ZWZ")
dt1_long = pivot_longer(dt1, New.method:ZWZ, 
                        names_to = "Method", values_to = 'value')
dt1_long$Method = factor(dt1_long$Method, levels=c("New.method", "Oracle", "ZWZ"))
p_Sn <- ggplot(data = dt1_long, aes(x=c1, y = value, colour = Method))+
  geom_line(aes(linetype=Method, size = Method))+theme_bw() +
  scale_color_manual(labels = c("New method", "Oracle", "Zhou et al's method"),values=c('black',"coral2", "cornflowerblue"))+
  scale_linetype_manual(labels = c("New method", "Oracle", "Zhou et al's method"),values=c("solid", "dashed", "twodash"))+
  theme(panel.grid =element_blank(),legend.position="bottom")+
  scale_size_manual(labels = c("New method", "Oracle", "Zhou et al's method"),values = c(0.7, 1.0,1.0))+
  scale_x_continuous(limits = c(-1.0,1.0),n.breaks = 11) + 
  labs( y="Empirical size and power", x= TeX("$c_1$")) + 
  geom_hline(yintercept = c(0.05), linetype = "dotted")
p_Sn

dt2 = df[22:42,c("c1","c2","R.Tn.Size.Power", "O.Tn.Size.Power","Z.Tn.Size.Power")]
colnames(dt2) = c("c1","c2", "New.method", "Oracle","ZWZ")
dt2_long = pivot_longer(dt2, New.method:ZWZ, names_to = "Method", values_to = 'value')
dt2_long$Method = factor(dt2_long$Method, levels=c("New.method", "Oracle", "ZWZ"))
p_Tn <- ggplot(data = dt2_long, aes(x=c2, y = value, colour = Method))+
  geom_line(aes(linetype=Method, size = Method))+theme_bw() +
  scale_color_manual(labels = c("New method", "Oracle", "Zhou et al's method"),values=c('black',"coral2", "cornflowerblue"))+
  scale_linetype_manual(labels = c("New method", "Oracle", "Zhou et al's method"), values=c("solid", "dashed","twodash"))+
  theme(panel.grid =element_blank())+
  scale_size_manual(labels = c("New method", "Oracle", "Zhou et al's method"), values = c(0.7, 1.0,1.0))+
  scale_x_continuous(limits = c(-1,1),n.breaks = 11) + 
  labs( y="Empirical size and power", x= TeX("$c_2$")) + 
  geom_hline(yintercept = c(0.05), linetype = "dotted")
p_Tn
ggarrange(p_Sn, p_Tn, ncol = 2, 
          common.legend = TRUE, legend = "bottom")
ggsave("results/confounder.pdf", dpi = 600, 
       width = 6.8, height =3,units = 'in')

## Generate latex table (to display the table correctly, 
## one should paste the output table to latex. Then add '&' in the heading
## and \begin{table} as well as \begin{tabular})


### Specify the rows to display
### Need to display c1 = -0.8,-0.4,0,0.4,0.8 when c2 = 0.5
###                 c2 = -0.8,-0.4,0,0.4,0.8 when c1 = 0.5
r_id = c(3,7,11,15,19, 24,28,32,36, 40) 
nr = length(r_id)
df = read.csv(summary_fn, header = T)
df3 = df[r_id,]

## Table 9: Estimated biases and standard deviations (in parentheses) of different methods with
## different c1 and c2 when confounding variables involved. Except for c1 and c2, the values in this
## table equal 100 times of the actual ones
Tab3 = data.frame(e1 = df3$c1,
                  e2 = rep('&', nr),
                  e3 = df3$c2,
                  e4 = rep('&', nr),
                  e5 = paste0('$',round(100*df3$R.alpha1.Bias, digits = 2),
                              '_{\\tiny{(',round(100*df3$R.alpha1.sample.std, digits = 2),')}}$'),
                  e6 = rep('&', nr),
                  e7 = paste0('$',round(100*df3$R.beta.Bias, digits = 2),
                              '_{\\tiny{(',round(100*df3$R.beta.sample.std, digits = 2),')}}$'),
                  e8 = rep('&', nr),
                  e9 = paste0('$',round(100*df3$O.alpha1.Bias, digits = 2),
                              '_{\\tiny{(',round(100*df3$O.alpha1.sample.std, digits = 2),')}}$'),
                  e10 = rep('&', nr),
                  e11 = paste0('$',round(100*df3$O.beta.Bias, digits = 2),
                               '_{\\tiny{(',round(100*df3$O.beta.sample.std, digits = 2),')}}$'),
                  e12 = rep('&', nr),
                  e13 = paste0('$',round(100*df3$Z.alpha1.Bias, digits = 2),
                                           '_{\\tiny{(',round(100*df3$Z.alpha1.sample.std, digits = 2),')}}$'),
                  e14 = rep('&', nr),
                  e15 = paste0('$',round(100*df3$Z.beta.Bias, digits = 2),
                                            '_{\\tiny{(',round(100*df3$Z.beta.sample.std, digits = 2),')}}$\\\\')
)
colnames(Tab3)<- c("$c_1$", " ", "$c_2$", " ", "$\\hat{\\alpha}_1$", " ", "$\\hat{\\beta}$",  " ", "$\\hat{\\alpha}^O_1$"," ", "$\\hat{\\beta}^O$",  " ","$\\hat{\\alpha}^Z_1$"," ","$\\hat{\\beta}^Z$")
write.table(Tab3, 'results/confounders_table9.txt', append = FALSE, sep = " ", dec = ".",
            row.names = FALSE, col.names = TRUE,quote =FALSE)

## Table 10: Estimated standard deviations and average estimated standard errors with their stan-
## dard deviations (in parentheses) over 500 replications with different c1 and c2 when confounding
## variables involved. Except for c1 and c2, the values in this table equal 100 times of the actual
## ones
Tab4 = data.frame(e1 = df3$c1,
                  e2 = rep('&', nr),
                  e3 = df3$c2,
                  e4 = rep('&', nr),
                  e5 = round(100*df3$R.alpha1.sample.std, digits = 2),
                  e6 = rep('&', nr),
                  e7 = paste0('$',round(100*df3$R.alpha1.est.std, digits = 2),
                              '_{\\tiny{(',round(100*df3$R.alpha.std.SE, digits = 2),')}}$'),
                  e8 = rep('&', nr),
                  e9 = round(100*df3$O.alpha1.sample.std, digits = 2),
                  e10 = rep('&', nr),
                  e11 = paste0('$',round(100*df3$O.alpha1.est.std, digits = 2),
                               '_{\\tiny{(',round(100*df3$O.alpha.std.SE, digits = 2),')}}$'),
                  e12 = rep('&', nr),
                  e13 = round(100*df3$R.beta.sample.std, digits = 2),
                  e14 = rep('&', nr),
                  e15 = paste0('$',round(100*df3$R.beta.est.std, digits = 2),
                               '_{\\tiny{(',round(100*df3$R.beta.std.SE, digits = 2),')}}$'),
                  e16 = rep('&', nr),
                  e17 = round(100*df3$O.beta.sample.std, digits = 2),
                  e18 = rep('&', nr),
                  e19 = paste0('$',round(100*df3$O.beta.est.std, digits = 2),
                               '_{\\tiny{(',round(100*df3$O.beta.std.SE, digits = 2),')}}$'),
                  e16 = rep('&', nr),
                  e17 = round(100*df3$Z.beta.sample.std, digits = 2),
                  e18 = rep('&', nr),
                  e19 = paste0('$',round(100*df3$Z.beta.est.std, digits = 2),
                               '_{\\tiny{(',round(100*df3$Z.beta.std.SE, digits = 2),')}}$\\\\')
                  
)
colnames(Tab4)<-c( "$c_1$", " ", "$c_2$", " ", "$\\hat{\\alpha}_1$ New method std", " ", "$\\hat{\\alpha}_1$ New method se(std)", " ", "$\\hat{\\alpha}_1$ Oracle std", " ", "$\\hat{\\alpha}_1$ Oracle se(std)",
                   " ", "$\\hat{\\beta}_1$ New method std", " ", "$\\hat{\\beta}_1$ New method se(std)", " ", "$\\hat{\\beta}_1$ Oracle std",  " ", "$\\hat{\\beta}_1$ Oracle se(std)", " ", "$\\hat{\\beta}_1$ Zhou std", " ", "$\\hat{\\beta}_1$ Zhou se(std)")
write.table(Tab4, 'results/confounder_table10.txt', append = FALSE, sep = " ", dec = ".",
            row.names = FALSE, col.names = TRUE,quote =FALSE)
