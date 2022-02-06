###########################################################################################
# Preprocessing for Section 4 simulation
load("results/preproc_JASA.Rdata")
source('utils.R')
## Standardize all variables
X = scale(X)
Y = scale(y)
M = scale(M)
S[,c(2,4:9)] = scale(S[,-c(1,3)])
S[,'sex'] = S[,'sex']-S[,'intercept']
S_imp = data.matrix(S[,c("intercept","sex")])

colnames(X) <- c("trauma")

## Measure the data dimension
n = length(X)
q = ncol(X)
p = ncol(M)
s = ncol(S)
# Tuning parameter lambda
lamb_min = 0.1
lamb_max = 0.228
ngrid = 50
lamb_grid = seq(lamb_min,lamb_max,length.out = ngrid)
hbic= c()
results =lapply(lamb_grid, HBIC_calc, xx=X, yy=Y, mm =M, S=S,n_imp = 8)
for( ii in 1: ngrid){
  hbic[ii] = results[[ii]]$BIC
}
plot(lamb_grid, hbic)
id = which(hbic==min(hbic))
id = tail(id,1)
alpha0_tld = results[[id]]$alpha0
lamb = lamb_grid[id]
print(paste("choose lambda:", lamb))
alpha0_hat = results[[id]]$alpha0
alpha1_hat = results[[id]]$alpha1
alpha2_hat = results[[id]]$alpha2
sigma1_hat = results[[id]]$sigma1_hat
## Selected mediators
A = which(alpha0_hat!=0)
alpha0_hat[A]
M_A = M[,A]
## Refit
W = data.matrix(cbind(X,S))
Gamma_hat = solve(t(W)%*%W) %*% t(W)%*% M
W_hat = W %*% Gamma_hat
W_hat_std = apply(W_hat, 2,sd)
res = M - W %*% Gamma_hat
res_std = apply(res, 2, sd)
refit = lm(Y~0+M_A+X+S)
alpha0_hat[A] = coef(refit)[1:length(A)]
alpha1_hat = coef(refit)['X']
alpha2_hat = coef(refit)[c('Sintercept','Sage','Ssex','Sbcell',
                           'SCD4T','SCD8T','Smonocytes','Sgranu','SNK' )]
dt = list(X = X,
          M = M,
          S = S,
          Y = Y)
coefs = list(alpha0 = alpha0_hat, alpha1=alpha1_hat,
             alpha2 = alpha2_hat, Gamma = Gamma_hat, sigma1 =summary(refit)$sigma)
save(dt, coefs, file = "results/simulation_allS.Rdata")

