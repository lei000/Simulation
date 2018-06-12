# Simulation
Envelope VAR model simulation

* library needed 
* code 
```library(MASS)

source(file="C:/Users/Lenovo/Desktop/Time series project/u.varenv.R")
source(file="C:/Users/Lenovo/Desktop/Time series project/u.varpenv.R")
source(file="C:/Users/Lenovo/Desktop/Time series project/varenv(2).R")
source(file="C:/Users/Lenovo/Desktop/Time series project/varpenv.R")
source(file="C:/Users/Lenovo/Desktop/Time series project/envlp/R/boot.env.R")
source(file="C:/Users/Lenovo/Desktop/Time series project/envlp/R/boot.penv.R")
source(file="C:/Users/Lenovo/Desktop/Time series project/envlp/R/boot.xenv.R")
source(file="C:/Users/Lenovo/Desktop/Time series project/envlp/R/contr.R")
source(file="C:/Users/Lenovo/Desktop/Time series project/envlp/R/cv.env.R")
source(file="C:/Users/Lenovo/Desktop/Time series project/envlp/R/cv.penv.R")
source(file="C:/Users/Lenovo/Desktop/Time series project/envlp/R/cv.xenv.R")
source(file="C:/Users/Lenovo/Desktop/Time series project/envlp/R/env.R")
source(file="C:/Users/Lenovo/Desktop/Time series project/envlp/R/envMU.R")
source(file="C:/Users/Lenovo/Desktop/Time series project/envlp/R/expan.R")
source(file="C:/Users/Lenovo/Desktop/Time series project/envlp/R/GE.R")
source(file="C:/Users/Lenovo/Desktop/Time series project/envlp/R/penv.R")
source(file="C:/Users/Lenovo/Desktop/Time series project/envlp/R/predict.env.R")
source(file="C:/Users/Lenovo/Desktop/Time series project/envlp/R/predict.penv.R")
source(file="C:/Users/Lenovo/Desktop/Time series project/envlp/R/predict.xenv.R")
source(file="C:/Users/Lenovo/Desktop/Time series project/envlp/R/predict2.env.R")
source(file="C:/Users/Lenovo/Desktop/Time series project/envlp/R/testcoef.env.R")
source(file="C:/Users/Lenovo/Desktop/Time series project/envlp/R/testcoef.penv.R")
source(file="C:/Users/Lenovo/Desktop/Time series project/envlp/R/testcoef.xenv.R")
source(file="C:/Users/Lenovo/Desktop/Time series project/envlp/R/u.env.R")
source(file="C:/Users/Lenovo/Desktop/Time series project/envlp/R/u.penv.R")
source(file="C:/Users/Lenovo/Desktop/Time series project/envlp/R/u.predict2.env.R")
source(file="C:/Users/Lenovo/Desktop/Time series project/envlp/R/u.xenv.R")
source(file="C:/Users/Lenovo/Desktop/Time series project/envlp/R/xenv.R")
source(file="C:/Users/Lenovo/Desktop/Time series project/envlp/R/data.trans.R")


N=100 # number of repetitions (simulation)
n = c(100, 200, 300, 500, 800, 1000, 1500) # sample size

K = length(n)
r = 5 # dimension of Y
m = 1 # dimension of Y
p = 5 # dimension of X
q = 1 # dimension of X, 1 in this case
u1 = 1  # envelope dimension u
sigma = .1
sigma0 = 1


#generate data
set.seed(123)
mu = rnorm(r)
eta1 = replicate(p, rnorm(u1))  
LL = replicate(r, runif(r))
LL = qr(LL) #qr square matrix
LL = qr.Q(LL, complete=TRUE)
L = LL[,1:u1]
L0 = LL[,(u1+1):r]


beta1 = L%*%t(eta1)
Sigma1 = L%*%t(L)*sigma^2+L0%*%t(L0)*sigma0^2
Sigma = Sigma1
Omega1 = sigma^2*diag(u1)
Omega10 = sigma0^2*diag(r-u1)




u_matrix = matrix(0, N, K)
var_pred_MSE = matrix(0, N, K)
en_pred_MSE = matrix(0, N, K)
beta_var=matrix(0, N, K)
beta_env=matrix(0, N, K)
ratio=matrix(0, N, K)
ratio_max=matrix(0, N, K)
ratio_min=matrix(0, N, K)
simu_coef=t(cbind(mu, beta1))

for (j in 1:K){
  
  for (k in 1:N){
    Y = matrix(0, r, n[j])
    Y[,1] = rnorm(r)
    epsil = t(mvrnorm(n[j], rep(0,r),Sigma))
    
    for (i in 2:n[j]) {
      Y[,i] = mu + beta1%*%Y[,i-1]+ epsil[,i]
    }
    
    data = t(Y)
    split_num = round((nrow(data))*0.8, 0)
    train = data[1: split_num,]
    test = data[(split_num+1):nrow(data) ,]
    
    library("MTS")
    #VARorder = VARorder(train)
    order = 1
    var_results <- MTS :: VAR(train, order)
    
    
    rownames(data) <- NULL
    Yts <- data[(order+1):nrow(data),]
    rownames(Yts) <- NULL
    Xts0 <- c(NA)
    for (i in 1:order){ 
      rownames(data[(order+1-i):(nrow(data)-i),]) <- NULL
      rownames(Xts0) <- NULL
      Xts0 <- cbind(Xts0,data[(order+1-i):(nrow(data)-i),]) 
      Xts <- Xts0[,c(-1)]
    }
    
    train_Xts <- Xts[1:(split_num-order), ]
    train_Yts <- Yts[1:(split_num-order), ]
    test_Xts <- Xts[(nrow(train_Xts)+1):nrow(Xts), ]
    test_Yts <- Yts[(nrow(train_Yts)+1):nrow(Yts), ]
    
    u_results = u.env(X=train_Xts, Y=train_Yts, alpha = 0.01)
    u_matrix[k, j] = u_results$u.bic
    varenv_results <- env(X=train_Xts, Y=train_Yts, u=u_results$u.bic, asy=TRUE)
    ratio[k, j]=mean(varenv_results$ratio)
    ratio_max[k, j]=max(varenv_results$ratio)
    ratio_min[k, j]=min(varenv_results$ratio)
    
    var_coef = var_results$coef
    #var_predict_value = rep(1, nrow(Xts)) %*% t(var_coef[1,]) + Xts %*% t(var_coef[2:nrow(var_coef),])
    var_predict_value = cbind(rep(1, nrow(test_Xts)), test_Xts)  %*% var_coef
    var_mean_square_error <- sqrt(mean ((var_predict_value - test_Yts)^2))
    var_pred_MSE[k, j] = var_mean_square_error
    beta_var[k, j]=sqrt(sum((var_coef-simu_coef)^2))
    #var_Y[, (k*4-3): (k*4)] = var_predict_value
    
    envar_coef = t(cbind(varenv_results$alpha, varenv_results$beta))
    #envar_predict_value = rep(1, nrow(Xts)) %*% t(envar_coef[1,]) + Xts%*% t(envar_coef[2:nrow(envar_coef),])
    envar_predict_value = cbind(rep(1, nrow(test_Xts)), test_Xts) %*% envar_coef
    envar_mean_square_error = sqrt(mean ((envar_predict_value - test_Yts)^2))
    en_pred_MSE[k, j] = envar_mean_square_error
    beta_env[k, j]=sqrt(sum((envar_coef-simu_coef)^2))

    
    
    
  }
}

bar_var_pred_MSE = apply(var_pred_MSE, 2, mean)
bar_env_pred_MSE = apply(en_pred_MSE, 2, mean)


bar_beta_var = apply(beta_var, 2, mean)
bar_beta_env = apply(beta_env, 2, mean)

bar_ratio=apply(ratio, 2, mean)
bar_ratio_max=apply(ratio_max, 2, max)
bar_ratio_min=apply(ratio_min, 2, min)

count_u <- function(x){
  u_count = sum(x==u1)
  return(paste(u_count,sep = "", "%"))
}
u_count = apply(u_matrix, 2, FUN = count_u)

#plot.ts(cbind(Yts[, 1], var_Y[, 1], en_Y[,1]), plot.type = "single", col=1:3)

plot(n, bar_var_pred_MSE, col="blue", pch=16, xlab="sample size", ylab="prediction error",
     ylim=c(min(bar_var_pred_MSE, bar_env_pred_MSE), max(bar_var_pred_MSE, bar_env_pred_MSE)))
lines(n, bar_var_pred_MSE, col="blue")
points(n, bar_env_pred_MSE, col="red", pch=16)
lines(n, bar_env_pred_MSE, col="red")
legend("topright",legend=c("VAR model","Envelope model"),
       text.col=c("blue","red"))

plot(n, bar_beta_var, col="blue", pch=16, xlab="sample size", ylab="estimated beta-true beta", 
     ylim=c(min(bar_beta_var,bar_beta_env), max(bar_beta_var,bar_beta_env)))
lines(n, bar_beta_var, col="blue")
points(n, bar_beta_env, col="red", pch=16)
lines(n, bar_beta_env, col="red")
legend("topright",legend=c("VAR model","Envelope model"),
       text.col=c("blue","red"))
```
