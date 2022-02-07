# Install or load package for Zhou, Wang \& Zhao (2020)
# Original freebird comes from https://github.com/rzhou14/freebird 
# but this package does not provide the estimated standard error of coefficients. 
# Instead, we modified its package output so that it can output standard error of the indirect effect beta
# Instead, we modified its package to output its standard error of the indirect effect beta
if(!require(freebird)){
  if(!require(devtools)){
    install.packages("devtools")
  }
  devtools::install_github("zengmudong/freebird")
  library(freebird)
}
if(!require(scalreg)){
  install.packages("scalreg")
}

if(!require(FDb.InfiniumMethylation.hg18)){
  if (!require(BiocManager)) {
    install.packages("BiocManager")
  }
  BiocManager::install("FDb.InfiniumMethylation.hg18")
  library(FDb.InfiniumMethylation.hg18)
}


library(glmnet)
options("digits" = 4)
LoadRawData<- function(InputPath, OutputPath = "./results/raw_data.Rdata"){
  #' Load raw data i tsv or txt format from specified path
  #' @param InputPath the absolute path of a folder that contain download data 
  #'                  from https://www.ebi.ac.uk/arrayexpress/experiments/E-GEOD-77445
  #' @param OutputPath the output path of merged data
  #' 
  #' @return OutputPath the output path of merged data          
  tab <- read.table(paste(InputPath,"X.tsv",sep =""),header = TRUE,sep="\t")
  # Extract x and y variables
  x <- tab$Characteristics..total.score.on.the.childhood.trauma.questionnaire.
  
  y <- tab$FactorValue..cortisol.stress.response.area.under.the.curve.auc.with.respect.to.the.increase.
  
  # Generate covariates: intercept, age, ...
  S <- data.frame(intercept = rep(1, nrow(tab)), 
                  age = tab$Characteristics..age., 
                  sex = tab$Characteristics..sex.,
                  bcell = tab$Characteristics..bcell.proportion.,
                  CD4T = tab$Characteristics..cd4.t.cell.proportion.,
                  CD8T = tab$Characteristics..cd8.t.cell.proportion.,
                  monocytes = tab$Characteristics..monocytes.cell.proportion.,
                  granu = tab$Characteristics..granulocytes.cell.proportion.,
                  NK = tab$Characteristics..natural.killer.cell.proportion.)
  S = data.matrix(S)
  
  # Load M variables
  indic <- sapply(strsplit(as.character(tab$Source.Name), " "), function(x) x[1])
  filename <- paste0(InputPath, indic[1], "_sample_table.txt")
  mdat <- read.delim(filename)
  M_raw <- matrix(0.0, nrow(tab), nrow(mdat))
  colnames(M_raw) <- mdat$Reporter.Identifier
  
  M_raw[1,] <- mdat$VALUE
  for (i in 2:length(indic)) { 
    filename <- paste0(InputPath, indic[i], "_sample_table.txt")
    M_raw[i,] <- read.delim(filename)$VALUE
  }
  print("Success of loading and merging data from .txt files")
  
  save(x, y, M_raw, S, file = OutputPath)
  
  print(paste0("Raw data is stored in: ",OutputPath))
  return(OutputPath)
}

MarginalScreening<-function(RawDataPath, topPotential = 1000){
  #' Marginal screening procedure based on Pearson correlation
  #' @param RawDataPath path of .Rdata that contain merged data
  #' @param topPotential number of top potential mediators will
  #'               be retained during screening step
  #'               
  #' @return OutputPath the path of .Rdata that stores 
  #'                    the marginal screening results 
  load(RawDataPath)
  ## Transform to M value; see Du et al., 2010
  M_log2 <- log2(M_raw / (1 - M_raw))
  
  
  ## 8 mediators suggested by Houtepen et al. (2016) and Kesteren & Oberski (2019)
  imp_m = c("cg27512205", "cg05608730", "cg26179948", "cg02309301",
            "cg12500973", "cg16657538", "cg25626453", "cg13136721")
  
  M_log2_ex = M_log2[,!colnames(M_log2) %in% imp_m]
  
  # Screening step to retainthe  top  1000  potential  mediators  
  # by  ranking  the  absolute  value  of  the  product  of  two  correlations  
  # -  between x and  each  element  of m,  and  
  # -  between y and  each  element  of m.   
  # This indeed is a marginal screening procedure based on Pearson correlation 
  # proposed by Fan & Lv(2008). 
  pcorr = abs(cor(M_log2_ex, y, method = 'pearson')
              *cor(M_log2_ex, x, method = 'pearson'))
  idxAB <- order(pcorr, decreasing = TRUE)[1:topPotential] 
  Mlog2_selAB <- M_log2_ex[, idxAB]
  # Combine the potential mediators with 8 mediators suggested 
  # by Houtepen et al. (2016) and Kesteren & Oberski (2019)
  M = cbind(Mlog2_selAB,M_log2[,imp_m])
  
  ## save result
  OutputPath = "./results/preproc_JASA.Rdata"
  X = as.matrix(x)
  colnames(X) <- "X"
  save(X, y, M, S, file = OutputPath)
  print(paste0("Preprocessed data is stored in ", OutputPath))
  return(OutputPath)
}

deSCAD <- function(z,lamb,a=3.7){
  # First order derivative of SCAD penalty
  # tuning parameter "a" (or "gamma") use the default value 3.7
  return(1*(z<=lamb)+pmax((a*lamb-z),0)/((a-1)*lamb)*(lamb<z))
}

## ONE-STEP SPARSE ESTIMATES IN NONCONCAVE PENALIZED LIKELIHOOD MODELS
LLA_h1<- function(X,Y,M,lamb,n,p,q, n_imp = 0,S = NULL){
  if(length(S) == 0){
    s = 0
    V = X
  }else{
    s = ncol(S)
    V = cbind(X,S)
  }
  # Step 1 using Lasso
  w1 = matrix(0,nrow = (p + q + s),ncol=1)
  w1[1:(p-n_imp)] = 1
  alpha_int = coef(glmnet(cbind(M,V),Y,family = 'gaussian', alpha=1,
                          lambda = lamb, penalty.factor=w1,intercept = FALSE))[-1]
  # Step 2 using linear approximation of SCAD
  w2 = matrix(0,nrow = (p+q + s),ncol=1)
  for(j in 1:(p- n_imp)){
    w2[j] = deSCAD(alpha_int[j],lamb)
  }
  alpha = coef(glmnet(cbind(M,V),Y,family = 'gaussian', alpha=1,
                      lambda = lamb, penalty.factor=w2, intercept = FALSE))
  return(alpha[-1])
}

HBIC_calc <- function(lamb, xx,yy,mm,S = NULL, n_imp =8){
  # Calculate HBIC for a specific tuning parameter lambda
  #' @param lamb a float value of tuning parameter lambda 
  #' @param xx The n by q exposure matrix. q can be 1, and q < n is required
  #' @param yy The n-dimensional outcome vector.
  #' @param mm The n by p mediator matrix. p can be larger than n.
  #' @param S The n by s confounding variables matrix. s can be 1, and s < n is required.
  #' @param n_imp an int for important mediators that will not be penalized
  
  #' @return
  #'        A list of class `HBIC_calc'
  #'        - BIC: HBIC score
  #'        - alpha0: estimated alpha0,
  #'        - alpha1: estimated alpha1,
  #'        - alpha2: estimated alpha2,
  #'        - sigma1_hat: estimated sigma1
  #'         
  n = nrow(xx)
  p = ncol(mm)
  q = ncol(xx)
  if(is.null(S)){
    s = 0
    result <- LLA_h1(xx,yy,mm,lamb,n,p,q, n_imp = n_imp)
    alpha0 = result[1:p]
    alpha1 = result[(p+1):(p+q)]
    alpha2 = NULL
    tmp = yy - mm%*%alpha0 - xx%*% alpha1
  }else{
    s = ncol(S)
    result <- LLA_h1(xx,yy,mm,lamb,n,p,q,n_imp = n_imp, S = S)
    alpha0 = result[1:p]
    alpha1 = result[(p+1):(p+q)]
    alpha2 = result[(p+q+1):(p+q+s)]
    
    tmp = yy - mm%*%alpha0 - xx%*% alpha1 - S %*% alpha2
  }
  
  df = length(which(alpha0!= 0))+q + s
  sigma_hat = t(tmp)%*%tmp/n
  BIC = log(sigma_hat) + df*log(log(n))*log(p+q + s)/n
  #obj = objective(xx,yy,M,alpha0,alpha1,lamb)
  return(list(BIC=BIC,alpha0=alpha0,alpha1 = alpha1, 
              alpha2 = alpha2, sigma1_hat = sigma_hat))
}


HBIC_PPLS <- function(X, y, M, S = NULL, lam_list =NULL, n_imp =8, topPotential= 1000){
  #' This function implements solving the partially penalized least squares with the SCAD penalty. 
  #' The tuning parameter lambda for the penalty function is chosen based on the high-dimensional 
  #' BIC (HBIC) method.
  #' @param X The n by q exposure matrix. q can be 1, and q < n is required
  #' @param y The n-dimensional outcome vector.
  #' @param M The n by p mediator matrix. p can be larger than n.
  #' @param S The n by s confounding variables matrix. s can be 1, and s < n is required. 
  #' @param lam_list a list of tuning parameter for HBIC
  #' @param topPotentialtop potential mediators that is retained from screening stage
  #' @return
  #'    A dataframe of selected important mediators (denoted as M_{\hat{A}} in the paper)
  hbic= c()
  
  result =lapply(lam_list, HBIC_calc, xx=X,yy=y, mm=M, S = S, n_imp = n_imp)
  
  for( ii in 1: length(lam_list)){
    hbic[ii] = result[[ii]]$BIC
  }
  # Find minimum HBIC score's corresponding lambda
  id = which(hbic==min(hbic))
  id = tail(id,1)
  lamb = lam_list[id]
  print(paste("choose lambda:", lamb))
  result = result[[id]]
  alpha0_hat = result$alpha0
  alpha1_hat = result$alpha1
  
  A = which(alpha0_hat!=0)
  # Selected mediators 
  M_A = M[,c(A[A>topPotential],A[A<=topPotential])]
  return(M_A)
}

MediatorExposure <- function(X,M){
  # Estimate \hat{Gamma} in M ~ Gamma X + Sigma
  
  Gamma_hat_sep = apply(M, 2, function(m){
    summary(lm(m~0+X))$coefficients[,'Estimate']
  })
  rownames(Gamma_hat_sep) <- colnames(X)
  Gamma_pvalue_sep = apply(M, 2, function(m){
    summary(lm(m~0+X))$coefficients[,'Pr(>|t|)']
  })
  Gamma_tvalue_sep = apply(M, 2, function(m){
    summary(lm(m~0+X))$coefficient[,'t value']
  })
  rownames(Gamma_pvalue_sep) <- colnames(X)
  return(list(Gamma_hat=t(Gamma_hat_sep),Gamma_pvalue=t(Gamma_pvalue_sep),
              Gamma_tvalue=t(Gamma_tvalue_sep)))
}

mediationInference<-function(X, Y, M, S){
  #' This function implements statistical inference of mediation model
  #' @param X The n by q exposure matrix. q can be 1, and q < n is required
  #' @param y The n-dimensional outcome vector.
  #' @param M The n by p selected important mediator matrix. p<n
  #' @param S The n by s confounding variables matrix. s can be 1, and s < n is required.
  #' 

  #' @return
  #'    A list of class `mediationInference`:
  #'    - beta_hat: estimated indirect effect
  #'    - alpha1_hat: estimated direct effect
  #'    - Sn: test statistc value of indirect effect
  #'    - Tn: test statistc value of direct effect
  #'    - summary_result: (dataframe) summary of The estimated coefficients, standard errors, test statistics values and p-values
  
  
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
    MS = cbind(M,S)
    alpha0_tld = solve(t(MS)%*%MS)%*%t(MS)%*%Y
    RSS02 = t(Y - MS%*% alpha0_tld) %*% (Y - MS%*%alpha0_tld)}
  
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
  Tn = (n-df) * (RSS02-RSS12)/RSS12
  
  p_beta = 1 - pchisq(Sn,1)
  p_alpha1 = 1 - pchisq(Tn,1)
  
  std_alpha1 = sqrt(var_alpha1_hat/n)
  std_beta =  sqrt(cov_beta_hat/n)
  options("digits" = 4)
  result.df = data.frame(Coeffcient = c("$\\mathbf{\\alpha}_1$", "$\\mathbf{\\beta}$"),
                         Estimated_Coeffcient = c(alpha1_hat, beta_hat),
                         SE =c(std_alpha1,std_beta),
                         Test_statistics = c(Tn, Sn),
                         p_value = c(p_alpha1, p_beta))
  
  return(list(Sn = Sn, Tn = Tn, 
              beta_hat = beta_hat, alpha0_hat = alpha0_hat,
              alpha1_hat = alpha1_hat, alpha2_hat = alpha2_hat, B = B, 
              var_beta = cov_beta_hat, var_alpha1_hat = var_alpha1_hat,
              p_beta = p_beta, p_alpha1 = p_alpha1,
              summary_result = result.df))
}


hdMediation <-function(X, y, M, lam_list, S = NULL, n_imp = 0){
  # This function implements the high-dimensional mediation analysis
  #' @param X (Required) n by q exposure matrix. q can be 1, and 
  #'                     q < n is required
  #' @param y (Required) n-dimensional outcome vector.
  #' @param M (Required) n by p mediator matrix. p can be larger than n.
  #' @param lam_list (Required) a list of tuning parameter for HBIC
  #' @param S (Optional) n by s confounding variables matrix. s can be 1, 
  #'                     and s < n is required. 
  #' @param n_imp (Optional) an int, specifying the number of unpenalized 
  #'                        mediators in the last n_imp columns of M matrix 
  #'                      
  #' @return
  #'    A list of class `hdMediation'
  #'    - Mediators_imp: Selected important mediators' name`
  #'    - summary_result: a data frame of estimated coefficients, 
  #'                      SE, Test_statistics	p_value
  
  M_A <- HBIC_PPLS(X, y, M, S = S, lam_list = lam_list, n_imp =n_imp)
  direct_indirect_eff <- mediationInference(X, y, M_A, S)
  
  
  return(list(Mediators_imp = colnames(M_A), 
              summary_result = direct_indirect_eff$summary_result))
}

analysis_helper <-function(X,Y,M, S = NULL, mod_name= NULL){
  
  ## This is a wrapper for mediation analysis with output of estimated 
  ## coefficients, test statistics, p-value and important mediators
  

  p = ncol(M)
  if(length(S) >0){s = ncol(S)}else{s = 0}
 
  A = c(1:ncol(M))
  M_A = M

  if(s ==0){
    df = data.frame(cbind(M_A, X))
    M_E = MediatorExposure(X,M_A)
  }
  else{
    V = cbind(X,S)
    df = data.frame(cbind(M_A, X, S))
    #Gamma_hat = solve(t(V)%*%V) %*%t(V)%*%M_A
    M_E = MediatorExposure(cbind(X,S),M_A)
  }
  fit4 = summary(lm(y~ 0+., data = df))
  
  Refit = mediationInference(X,Y,M_A,S)
  
  mod_summary = data.frame(round(fit4$coefficients[,c(1,2,4)],digits = 3))
  if(!is.null(mod_name)){
    colnames(mod_summary) <- paste(mod_name, c("Estimate", "SE", "p-value"))
  }else{colnames(mod_summary) <- c("Estimate", "SE", "p-value")}
  

  return(list(mod_summary = mod_summary, 
              summary_result = Refit$summary_result, A_hat = A,
              Gamma_hat= round(M_E$Gamma_hat,digits = 3),
              Gamma_pvalue= round(M_E$Gamma_pvalue, digits = 3),
              Gamma_tvalue= round(M_E$Gamma_tvalue,digits = 3)
  ))
}


Fisher.corr.test<- function(r, n){
  if(r == 1){fisherZ = Inf
  return(0)
  }
  else{
    fisherZ = abs(0.5*log((1+r)/(1-r)))
    return(2*(1 - pnorm(fisherZ, 0, 1/sqrt(n-3))))
  }
}

Fisher.helper<-function(df){
  #' Calculate Pearson correlation and p-value of Fisher correlation test
  #' @param df dataframe n * p containing the variables for Fisher test
  #' @return output dataframe, lower triangle is Fisher correlation, 
  #'                the upper triagnle is p-value of Fisher test
  n = nrow(df)
  cor_mat = cor(df)
  p.value = apply(cor_mat, c(1,2), Fisher.corr.test, n =n)
  
  output <- cor_mat
  output[upper.tri(output)] <- p.value[upper.tri(p.value)]
  output = format(round(output, 4), nsmall = 4)
  diag(output) <- ''
  return(output)
}

partial.helper <-function(df, xy, zz){
  #' Calculate marginal corr between X and Y given Z
  #' @param df dataframe n * p containing the variables for Fisher test
  #' @param xy list column id for correlation test
  #' @param zz list column id for variables to be conditioned
  #' @return output dataframe, lower triangle is Fisher correlation, 
  #'                the upper triagnle is p-value of Fisher test
  p = length(xy)
  estimate = matrix(rep(1,p*p),p,p)
  colnames(estimate) = colnames(df)[xy]
  rownames(estimate) = colnames(df)[xy]
  p.value = matrix(rep(0,p*p),p,p)
  colnames(p.value) = colnames(df)[xy]
  rownames(p.value) = colnames(df)[xy]
  for(x.id in 1:(p - 1)){
    for(y.id in (x.id+1):p){
      pcor_pearson = ppcor::pcor(df[,c(xy[x.id], xy[y.id], zz)], method = "pearson")
      estimate[x.id, y.id] = pcor_pearson$estimate[1,2]
      estimate[y.id, x.id] = pcor_pearson$estimate[1,2]
      p.value[x.id, y.id] = pcor_pearson$p.value[1,2]
      p.value[y.id, x.id] = pcor_pearson$p.value[1,2]
    }
  }
  
  output <- estimate
  output[upper.tri(output)] <- p.value[upper.tri(p.value)]
  output = format(round(output, 4), nsmall = 4)
  diag(output) <- ''
  return(output)
}

generateTableS4<-function(models){
  R2 = c()
  F_stat = c()
  p_val = c()
  for(mod in models){
    R2 = append(R2, mod$r.squared)
    F_stat = append(F_stat, mod$fstatistic[1])
    p_val = append(p_val, pf(mod$fstatistic[1],mod$fstatistic[2],mod$fstatistic[3],lower.tail=FALSE))
  }
  return(list(R2 = R2, F_stat = F_stat, p_val = p_val))
}

