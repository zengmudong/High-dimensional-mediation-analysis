Real Data Analysis
================
Xu Guo, Runze Li, Jingyuan Liu, Mudong Zeng

This Rmd is devoted to an empirical analysis of the same data set as
that in Houtepen et al. (2016) and van Kesteren & Oberski (2019), for
studying how DNA methylation plays a role in the regulation of human
stress reactivity. This Rmd aims to reproduce the results in Section 3.

The data can be downloaded from the following website:
<https://www.ebi.ac.uk/arrayexpress/experiments/E-GEOD-77445>, and the
data set consists of 385,882 DNA methylation loci and various variables
for 85 people.

# Preprocess raw data

The original data available to the authors are .txt files. We loaded
these files and merged into a clean data set (containing two vectors,
exposure variable X as well as response variable y, and two dataframes
potential mediators M as well as confounding variables S). To make it
convenient for further analysis, we stored it in a .RData file.

## Load data frame from participant info

## Marginal screening for potential mediators

We carry out a screening step to retain the top 1000 potential mediators
by ranking the absolute value of the product of two correlations - the
correlation between
![\\mathbf{x}](https://latex.codecogs.com/png.latex?%5Cmathbf%7Bx%7D
"\\mathbf{x}") and each element of
![\\mathbf{m}](https://latex.codecogs.com/png.latex?%5Cmathbf%7Bm%7D
"\\mathbf{m}"), and between ![y](https://latex.codecogs.com/png.latex?y
"y") and each element of
![\\mathbf{m}](https://latex.codecogs.com/png.latex?%5Cmathbf%7Bm%7D
"\\mathbf{m}"). This indeed is a marginal screening procedure based on
Pearson correlation proposed by Fan & Lv (2008).

``` r
RawDataPath = "results/raw_data.Rdata"
topPotential = 1000
cleanedData = MarginalScreening(RawDataPath, topPotential = topPotential)
```

    ## [1] "Preprocessed data is stored in ./results/preproc_JASA.Rdata"

``` r
load(cleanedData)
```

# Our proposed new method

Procedure of identifying active mediators in the mediation models (2.1)
and (2.2), and estimating the direct effect
![\\mathbf{\\alpha\_1}](https://latex.codecogs.com/png.latex?%5Cmathbf%7B%5Calpha_1%7D
"\\mathbf{\\alpha_1}") and indirect effect
![\\boldsymbol{\\beta}](https://latex.codecogs.com/png.latex?%5Cboldsymbol%7B%5Cbeta%7D
"\\boldsymbol{\\beta}") that can get around high dimensional matrix
estimation.

## Identifying active mediators

we apply the partial penalized least squared method to fit model (2.1)
by only penalizing
![\\boldsymbol{\\alpha}\_0](https://latex.codecogs.com/png.latex?%5Cboldsymbol%7B%5Calpha%7D_0
"\\boldsymbol{\\alpha}_0"), which corresponding to Eq. (2.5).

### Algorithm for solving Eq. (2.5)

We apply the local linear approximation algorithm (LLA) in Zou and Li
(2008) with the SCAD penalty (Fan & Li, 2001) ,   
![&#10;p'\_{\\lambda}(t)=\\lambda\\{I(t\\leq
\\lambda)&#10;+\\frac{(a\\lambda-t)\_{+}}{(a-1)\\lambda}I(t\>\\lambda)\\},&#10;](https://latex.codecogs.com/png.latex?%0Ap%27_%7B%5Clambda%7D%28t%29%3D%5Clambda%5C%7BI%28t%5Cleq%20%5Clambda%29%0A%2B%5Cfrac%7B%28a%5Clambda-t%29_%7B%2B%7D%7D%7B%28a-1%29%5Clambda%7DI%28t%3E%5Clambda%29%5C%7D%2C%0A
"
p'_{\\lambda}(t)=\\lambda\\{I(t\\leq \\lambda)
+\\frac{(a\\lambda-t)_{+}}{(a-1)\\lambda}I(t\>\\lambda)\\},
")  
and set ![a=3.7](https://latex.codecogs.com/png.latex?a%3D3.7 "a=3.7").
The tuning parameter
![\\lambda](https://latex.codecogs.com/png.latex?%5Clambda "\\lambda")
for our method is chosen based on the high-dimensional BIC (HBIC) method
in Wang et al. (2013). For a fixed regularization parameter
![\\lambda](https://latex.codecogs.com/png.latex?%5Clambda "\\lambda"),
define   
![(\\hat{\\boldsymbol{\\alpha}}\_0^{\\lambda},
\\hat{\\boldsymbol{\\alpha}}\_1^{\\lambda})=\\min\_{\\boldsymbol{\\alpha}\_0,\\boldsymbol{\\alpha}\_1}\\frac{1}{2n}\\|\\mathbf{Y}-\\mathbf{M}\\boldsymbol{\\alpha}\_{0}-\\mathbf{X}\\boldsymbol{\\alpha}\_1\\|^2\_2+\\sum\_{j=1}^p
p\_{\\lambda}(|\\alpha\_{0,j}|).](https://latex.codecogs.com/png.latex?%28%5Chat%7B%5Cboldsymbol%7B%5Calpha%7D%7D_0%5E%7B%5Clambda%7D%2C%20%5Chat%7B%5Cboldsymbol%7B%5Calpha%7D%7D_1%5E%7B%5Clambda%7D%29%3D%5Cmin_%7B%5Cboldsymbol%7B%5Calpha%7D_0%2C%5Cboldsymbol%7B%5Calpha%7D_1%7D%5Cfrac%7B1%7D%7B2n%7D%5C%7C%5Cmathbf%7BY%7D-%5Cmathbf%7BM%7D%5Cboldsymbol%7B%5Calpha%7D_%7B0%7D-%5Cmathbf%7BX%7D%5Cboldsymbol%7B%5Calpha%7D_1%5C%7C%5E2_2%2B%5Csum_%7Bj%3D1%7D%5Ep%20p_%7B%5Clambda%7D%28%7C%5Calpha_%7B0%2Cj%7D%7C%29.
"(\\hat{\\boldsymbol{\\alpha}}_0^{\\lambda}, \\hat{\\boldsymbol{\\alpha}}_1^{\\lambda})=\\min_{\\boldsymbol{\\alpha}_0,\\boldsymbol{\\alpha}_1}\\frac{1}{2n}\\|\\mathbf{Y}-\\mathbf{M}\\boldsymbol{\\alpha}_{0}-\\mathbf{X}\\boldsymbol{\\alpha}_1\\|^2_2+\\sum_{j=1}^p p_{\\lambda}(|\\alpha_{0,j}|).")  
The minimization of the partially penalized least squares method can be
carried out as follows.

1.  Get initial values for
    ![\\boldsymbol{\\alpha}^{(0)}\_0,\\boldsymbol{\\alpha}^{(0)}\_1](https://latex.codecogs.com/png.latex?%5Cboldsymbol%7B%5Calpha%7D%5E%7B%280%29%7D_0%2C%5Cboldsymbol%7B%5Calpha%7D%5E%7B%280%29%7D_1
    "\\boldsymbol{\\alpha}^{(0)}_0,\\boldsymbol{\\alpha}^{(0)}_1") by
    minimizing a partial
    ![L\_1](https://latex.codecogs.com/png.latex?L_1 "L_1")-penalized
    least squares:   
    ![(\\hat{\\boldsymbol{\\alpha}}\_0^{(0)},
    \\hat{\\boldsymbol{\\alpha}}\_1^{(0)})=\\min\_{\\boldsymbol{\\alpha}\_0,\\boldsymbol{\\alpha}\_1}&#10;\\frac{1}{2n}\\|\\mathbf{Y}-\\mathbf{M}\\boldsymbol{\\alpha}\_{0}-\\mathbf{X}\\boldsymbol{\\alpha}\_1\\|^2\_2+\\lambda\\sum\_{j=1}^p
    |\\alpha\_{0,j}|.](https://latex.codecogs.com/png.latex?%28%5Chat%7B%5Cboldsymbol%7B%5Calpha%7D%7D_0%5E%7B%280%29%7D%2C%20%5Chat%7B%5Cboldsymbol%7B%5Calpha%7D%7D_1%5E%7B%280%29%7D%29%3D%5Cmin_%7B%5Cboldsymbol%7B%5Calpha%7D_0%2C%5Cboldsymbol%7B%5Calpha%7D_1%7D%0A%5Cfrac%7B1%7D%7B2n%7D%5C%7C%5Cmathbf%7BY%7D-%5Cmathbf%7BM%7D%5Cboldsymbol%7B%5Calpha%7D_%7B0%7D-%5Cmathbf%7BX%7D%5Cboldsymbol%7B%5Calpha%7D_1%5C%7C%5E2_2%2B%5Clambda%5Csum_%7Bj%3D1%7D%5Ep%20%7C%5Calpha_%7B0%2Cj%7D%7C.
    "(\\hat{\\boldsymbol{\\alpha}}_0^{(0)}, \\hat{\\boldsymbol{\\alpha}}_1^{(0)})=\\min_{\\boldsymbol{\\alpha}_0,\\boldsymbol{\\alpha}_1}
\\frac{1}{2n}\\|\\mathbf{Y}-\\mathbf{M}\\boldsymbol{\\alpha}_{0}-\\mathbf{X}\\boldsymbol{\\alpha}_1\\|^2_2+\\lambda\\sum_{j=1}^p |\\alpha_{0,j}|.")  

2.  Solve   
    ![(\\hat{\\boldsymbol{\\alpha}}\_0^{(k+1)},
    \\hat{\\boldsymbol{\\alpha}}\_1^{(k+1)})=\\min\_{\\boldsymbol{\\alpha}\_0,\\boldsymbol{\\alpha}\_1}&#10;\\frac{1}{2n}\\|\\mathbf{Y}-\\mathbf{M}\\boldsymbol{\\alpha}\_{0}-\\mathbf{X}\\boldsymbol{\\alpha}\_1\\|^2\_2+\\sum\_{j=1}^p
    p'\_{\\lambda}(|\\alpha^{(k)}\_{0,j}|)|\\alpha\_{0,j}|, \\text{ for
    }
    k=1,2,\\cdots,](https://latex.codecogs.com/png.latex?%28%5Chat%7B%5Cboldsymbol%7B%5Calpha%7D%7D_0%5E%7B%28k%2B1%29%7D%2C%20%5Chat%7B%5Cboldsymbol%7B%5Calpha%7D%7D_1%5E%7B%28k%2B1%29%7D%29%3D%5Cmin_%7B%5Cboldsymbol%7B%5Calpha%7D_0%2C%5Cboldsymbol%7B%5Calpha%7D_1%7D%0A%5Cfrac%7B1%7D%7B2n%7D%5C%7C%5Cmathbf%7BY%7D-%5Cmathbf%7BM%7D%5Cboldsymbol%7B%5Calpha%7D_%7B0%7D-%5Cmathbf%7BX%7D%5Cboldsymbol%7B%5Calpha%7D_1%5C%7C%5E2_2%2B%5Csum_%7Bj%3D1%7D%5Ep%20p%27_%7B%5Clambda%7D%28%7C%5Calpha%5E%7B%28k%29%7D_%7B0%2Cj%7D%7C%29%7C%5Calpha_%7B0%2Cj%7D%7C%2C%20%5Ctext%7B%20for%20%7D%20k%3D1%2C2%2C%5Ccdots%2C
    "(\\hat{\\boldsymbol{\\alpha}}_0^{(k+1)}, \\hat{\\boldsymbol{\\alpha}}_1^{(k+1)})=\\min_{\\boldsymbol{\\alpha}_0,\\boldsymbol{\\alpha}_1}
\\frac{1}{2n}\\|\\mathbf{Y}-\\mathbf{M}\\boldsymbol{\\alpha}_{0}-\\mathbf{X}\\boldsymbol{\\alpha}_1\\|^2_2+\\sum_{j=1}^p p'_{\\lambda}(|\\alpha^{(k)}_{0,j}|)|\\alpha_{0,j}|, \\text{ for } k=1,2,\\cdots,")  
        until ![\\{(\\hat{\\boldsymbol{\\alpha}}\_0^{(k)},
    \\hat{\\boldsymbol{\\alpha}}\_1^{(k)})\\}](https://latex.codecogs.com/png.latex?%5C%7B%28%5Chat%7B%5Cboldsymbol%7B%5Calpha%7D%7D_0%5E%7B%28k%29%7D%2C%20%5Chat%7B%5Cboldsymbol%7B%5Calpha%7D%7D_1%5E%7B%28k%29%7D%29%5C%7D
    "\\{(\\hat{\\boldsymbol{\\alpha}}_0^{(k)}, \\hat{\\boldsymbol{\\alpha}}_1^{(k)})\\}")
    converges.

Below is the implementation of solving Eq. (2.5). The candidate set of
![\\lambda](https://latex.codecogs.com/png.latex?%5Clambda "\\lambda")
is 50 equally spaced grid from 20 to 70.

``` r
load("results/preproc_JASA.Rdata")
topPotential = 1000

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
  # S is cofounder
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


HBIC_PPLS <- function(X, y, M, S = NULL, lam_list =NULL, n_imp =8){
    #' This function implements solving the partially penalized least squares with the SCAD penalty. The tuning parameter lambda for the penalty function is chosen based on the high-dimensional BIC (HBIC) method.
  #' @param X The n by q exposure matrix. q can be 1, and q < n is required
  #' @param y The n-dimensional outcome vector.
  #' @param M The n by p mediator matrix. p can be larger than n.
  #' @param S The n by s confounding variables matrix. s can be 1, and s < n is required.
  #' 
  #' @param lam_list a list of tuning parameter for HBIC
  
  #' @return
  #'
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
  
}

# specify the candidate set of lambda
lamb_min = 20
lamb_max = 70
ngrid = 50
lam_list = seq(lamb_min,lamb_max,length.out = ngrid)

M_A = HBIC_PPLS(X, y, M, S = S, lam_list =lam_list, n_imp =8)
```

    ## [1] "choose lambda: 60.8163265306122"

### Estimated Coefficients, SE, t-values and p-values (Table 2)

``` r
# Refit to produce Table 2:
colnames(S) = paste0("Z",0:8)
df = data.frame(cbind(M_A, X, S))
round(summary(lm(y~ 0+., data = df))$coefficients,digits = 3)
```

    ##             Estimate Std. Error t value Pr(>|t|)
    ## cg27512205  -237.547    199.506  -1.191    0.238
    ## cg05608730  -301.168    151.038  -1.994    0.050
    ## cg26179948  -474.486    160.042  -2.965    0.004
    ## cg02309301   259.730    108.633   2.391    0.020
    ## cg12500973    30.029    116.354   0.258    0.797
    ## cg16657538    84.330     53.236   1.584    0.118
    ## cg25626453   369.183     97.988   3.768    0.000
    ## cg13136721   260.990     65.585   3.979    0.000
    ## cg19230917   321.196    149.918   2.142    0.036
    ## cg06422529   418.173    107.252   3.899    0.000
    ## cg03199124   471.865    143.943   3.278    0.002
    ## X              1.365      4.553   0.300    0.765
    ## Z0         -3110.834   3805.517  -0.817    0.417
    ## Z1            -1.864      2.056  -0.906    0.368
    ## Z2           349.037     82.811   4.215    0.000
    ## Z3          1843.451   3702.100   0.498    0.620
    ## Z4           406.642   3533.801   0.115    0.909
    ## Z5           781.938   3368.749   0.232    0.817
    ## Z6           967.962   3745.283   0.258    0.797
    ## Z7           123.714   3544.597   0.035    0.972
    ## Z8           341.974   3401.318   0.101    0.920

### Estimated Coefficients of ![\\Gamma\_1](https://latex.codecogs.com/png.latex?%5CGamma_1 "\\Gamma_1") and ![\\Gamma\_2](https://latex.codecogs.com/png.latex?%5CGamma_2 "\\Gamma_2") and their p-values (Table 3)

``` r
colnames(M_A) = paste0("m",1:11)
Gamma = MediatorExposure(cbind(X, S), M_A)
round(Gamma$Gamma_hat,digits = 3)
```

    ##          X     Z0     Z1     Z2     Z3     Z4     Z5     Z6     Z7     Z8
    ## m1   0.005 -2.999  0.000  0.016  0.870  0.197 -0.311  0.784  0.406  0.010
    ## m2   0.007 -1.723 -0.002  0.013  1.479  0.032  0.693  0.602  0.935  1.448
    ## m3   0.006 -3.566 -0.001  0.004  0.985  0.704 -0.809  0.840  0.567 -0.373
    ## m4  -0.012 -9.222  0.007 -0.115  5.179  3.716  4.666  5.653  4.371  4.901
    ## m5  -0.011 -3.974  0.004 -0.009  2.051 -2.022 -0.485  0.437 -0.860 -0.594
    ## m6   0.023  0.322  0.003  0.403 -6.719 -5.542  0.115 -1.608 -5.701 -4.943
    ## m7  -0.012  5.701  0.002  0.000  0.562 -0.290  0.552 -0.510 -0.285 -0.516
    ## m8   0.012  1.093  0.001  0.083 -2.086  3.799  1.748  3.922  2.673  1.312
    ## m9  -0.006 -2.009  0.001  0.018  1.257 -1.623 -1.448 -0.942 -1.181 -1.276
    ## m10 -0.008  8.459  0.003  0.042 -7.930 -4.387 -2.036 -3.341 -4.408 -4.387
    ## m11 -0.006 -5.736  0.003 -0.080  4.745  2.151  1.013  1.080  3.009  2.688

``` r
round(Gamma$Gamma_pvalue,digits = 3)
```

    ##         X    Z0    Z1    Z2    Z3    Z4    Z5    Z6    Z7    Z8
    ## m1  0.044 0.225 0.711 0.766 0.726 0.936 0.895 0.765 0.868 0.997
    ## m2  0.016 0.568 0.145 0.848 0.627 0.992 0.811 0.851 0.755 0.613
    ## m3  0.020 0.201 0.448 0.946 0.724 0.798 0.761 0.776 0.837 0.887
    ## m4  0.004 0.025 0.001 0.202 0.205 0.356 0.230 0.191 0.277 0.203
    ## m5  0.002 0.286 0.054 0.917 0.584 0.584 0.892 0.912 0.816 0.866
    ## m6  0.004 0.968 0.490 0.026 0.406 0.487 0.988 0.850 0.474 0.515
    ## m7  0.006 0.205 0.444 0.998 0.901 0.948 0.898 0.915 0.949 0.903
    ## m8  0.056 0.865 0.861 0.563 0.748 0.554 0.777 0.568 0.676 0.830
    ## m9  0.041 0.528 0.551 0.795 0.695 0.608 0.635 0.781 0.709 0.672
    ## m10 0.048 0.044 0.223 0.646 0.060 0.288 0.608 0.449 0.286 0.266
    ## m11 0.045 0.068 0.108 0.253 0.133 0.488 0.734 0.744 0.332 0.364

## Test of direct effect and indirect effect

### The estimated coefficients, standard errors, test statistics values and p-values (Table 4)

``` r
direct_indirect_eff <- mediationInference(X, y, M_A, S)
direct_indirect_eff$summary_result
```

    ##    Coeff Estimated_Coeffcient   std Test_statistics  p_value
    ## 1 alpha1                1.365 4.553         0.08993 0.764268
    ## 2   beta              -17.373 5.494         9.99708 0.001568
