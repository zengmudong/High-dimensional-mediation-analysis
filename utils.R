# Install or load package from Zhou, Wang \& Zhao (2020)
# Origial freebird in comes from https://github.com/rzhou14/freebird 
# but this package does not provide the estimated standard error of coefficients. 
# Instead, we modified its package output so that it can output standard error of the indirect effect beta

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

library(scalreg)
library(glmnet)
options("digits" = 4)
LoadRawData<- function(InputPath, OutputPath = "./results/raw_data.Rdata"){
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
  # topPotential: number of top potential mediators will 
  #               be retained during screening step
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
  result.df = data.frame(Coeff = c("alpha1", "beta"),
                         Estimated_Coeffcient = c(alpha1_hat, beta_hat),
                         std =c(std_alpha1,std_beta),
                         Test_statistics = c(Tn, Sn),
                         p_value = c(p_alpha1, p_beta))
  
  return(list(Sn = Sn, Tn = Tn, 
              beta_hat = beta_hat, alpha0_hat = alpha0_hat,
              alpha1_hat = alpha1_hat, alpha2_hat = alpha2_hat, B = B, 
              var_beta = cov_beta_hat, var_alpha1_hat = var_alpha1_hat,
              p_beta = p_beta, p_alpha1 = p_alpha1,
              summary_result = result.df))
}
