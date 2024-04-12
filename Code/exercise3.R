lib <- c('fBasics','data.table', 'xts','zoo', 'quantmod', 'PerformanceAnalytics', 
         'PortfolioAnalytics', 'GenSA', 'timeSeries','fPortfolio', 'quantmod', 
         'plyr', 'PortfolioAnalytics', 'ROI', 'GenSA','DEoptim', 'xtable', 
         "ROI.plugin.glpk", 'ROI.plugin.quadprog','portfolioBacktest','ggplot2','mvtnorm')
loading.lib <- lapply(lib, require, character.only = TRUE)
#library(mvtnorm)
## 1-	Select 9 stocks from the dataset:
fileloc <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(fileloc)
rm(fileloc)
dataset <- readRDS("data/dataset.rds")
dataset$PNlog

## 2-	Retrieve Factor Models:
# upload Fama-French factors
mydata <- read.csv("data/F-F_Research_Data_Factors_daily.CSV", skip = 4)
mydata <- mydata[-nrow(mydata), ]  # remove last row
fama_lib <- xts(x = mydata[, c(2,3,4)], order.by = as.Date(paste(mydata[, 1]), "%Y%m%d"))
str(fama_lib)
fama_lib

begin_date <- "2018-01-01"   #"2015-01-01"
end_date <- "2020/12/20"     #"2019-12-31" #"2017-12-31"  #
period<-paste(begin_date,"/",end_date,sep="")

## 3-	Calculate Factor Returns:
# prepare stock data
#colnames(dataset$adjusted)
stockPrices<-dataset$adjusted[period]     #[1:1250]
stockPrices <- stockPrices[, 1:9]
## for stocks "AAPL" "ABBV" "AMZN" "DB"   "DIS"  "FB" "GOOG" "HAL" "HSBC" - 9 stocks
tclass(stockPrices) <- "Date"
X <- diff(log(stockPrices), na.pad = FALSE)
N <- ncol(X)  # number of stocks
T <- nrow(X)  # number of days

# prepare Fama-French factors
F_FamaFrench <- fama_lib[index(X)]/100
F_FamaFrench

##create the PNlog Market index
PNlogMkt <- na.approx(as.xts(rowAvgs(dataset$PNlog),order.by=index(dataset$PNlog)))


SentIndx <- PNlogMkt[index(X)]

## 4-	Construct the Covariance Matrix (Sigma):
# sample covariance matrix
#Sigma_SCM <- cov(X) #?????


### First model
# Fama-French 3-factor model
F_ <- cbind(ones = 1, F_FamaFrench)
Gamma <- t(solve(t(F_) %*% F_, t(F_) %*% X)) 

# Gamma - contiene los coeficientes estimados del modelo de tres factores de Fama-French
colnames(Gamma) <- c("alpha", "beta1", "beta2", "beta3")
alpha <- Gamma[, 1]
B <- Gamma[, 2:4]

# Combine them into a matrix for all factors
FF_returns <- cbind(alpha, B)
# Construct the Sigma Matrix using FF_returns:
Sigma_3FF <- cov(FF_returns)

E <- xts(t(t(X) - Gamma %*% t(F_)), index(X)) # Compute Residuals
PsiFF <- (1/(T-4)) * t(E) %*% E # Compute Residual Covariance
Sigma_FamaFrench <- B %*% cov(F_FamaFrench) %*% t(B) + diag(diag(PsiFF)) # needed for computing the Robust Global Maximum Return Portfolio
Sigma_FamaFrench

### Second model
# 1-factor model SentInx
F_ <- cbind(ones = 1, SentIndx)
Gamma <- t(solve(t(F_) %*% F_, t(F_) %*% X))
colnames(Gamma) <- c("alpha", "beta")
alpha <- Gamma[, 1]
beta <- Gamma[, 2]
E <- xts(t(t(X) - Gamma %*% t(F_)), index(X))
Psi_Sent <- (1/(T-2)) * t(E) %*% E
Sigma_SentInx <- as.numeric(var(SentIndx)) * beta %o% beta + diag(diag(Psi_Sent))

## 5-	Compute the Robust (Ellipsoid) Global Maximum Return Portfolio:
##############################################
## Robust optimization of global maximum return portfolio
## Problem: max_w: mu*w ; subject to: w1+w2+...+wn=1, w>=0
## We have seen optimal solution almost always allocate all the budget 
## on the single stock with largest expected return. 
## To prevent this assume mu belongs to some convex uncertainty set:
## mu' = colMeans(X) + k*Sigma^{0.5}*u ,  ||u||< 1, and find minimum such w*mu'
## This minimum is w*mu - k*||Sigma^{0.5}*w|| 
## Implement solution with a solver (CVXR)

library(CVXR)

mu <- colMeans(X)
portfolioMaxReturnRobustEllipsoid <- function(mu_hat, S, kappa = 0.1) {
  S12 <- chol(S)  # t(S12) %*% S12 = Sigma
  w <- Variable(length(mu_hat))
  prob <- Problem(Maximize( t(w) %*% mu_hat - kappa*p_norm(S12 %*% w,p=2) ), 
                  constraints = list(w >= 0, sum(w) == 1))
  result <- solve(prob)
  return(as.vector(result$getValue(w)))
}

# A robust solution
kappa<-0.9
w_GMRP_robust <- portfolioMaxReturnRobustEllipsoid(mu,Sigma_FamaFrench,kappa)
names(w_GMRP_robust) <- colnames(X)
w_all_GMRP_robust_ellipsoid <- cbind(w_GMRP_robust)

dev.off()
# plot allocation
barplot(t(w_GMRP_robust), col = "red", legend = colnames(X), beside = TRUE,
        main = paste("Robust (ellipsoid) Global maximum return portfolio allocation, k=",kappa),
        xlab = "stocks", ylab = "dollars")


# Try multiple robust noisy solutions to check for sensitivity.
set.seed(357)
for (i in 1:6) {
  X_noisy <- rmvnorm(n = T, mean = mu, sigma = Sigma_FamaFrench)
  mu_noisy <- colMeans(X_noisy)
  Sigma_noisy <- cov(X_noisy)
  
  w_GMRP_robust_ellipsoid_noisy <- portfolioMaxReturnRobustEllipsoid(mu_noisy, Sigma_noisy, kappa)
  w_all_GMRP_robust_ellipsoid <- cbind(w_all_GMRP_robust_ellipsoid, w_GMRP_robust_ellipsoid_noisy)
}

# plot to compare the allocations
barplot(t(w_all_GMRP_robust_ellipsoid), col = 1:7, #legend = colnames(w_all_GMRP_robust_ellipsoid), 
        beside = TRUE, args.legend = list(bg = "white"),
        main = paste("Robust (ellipsoid) Global maximum return portfolio allocation, k=",kappa), xlab = "stocks", ylab = "dollars")


chart.StackedBar(t(w_all_GMRP_robust_ellipsoid), 
                 main = paste("Robust (ellipsoid) Global maximum return portfolio allocation, k=",kappa), 
                 ylab = "w", space = 0, border = NA)

####################
####################
####################

# A robust solution
kappa<-0.63
w_GMRP_robust <- portfolioMaxReturnRobustEllipsoid(mu,Sigma_SentInx,kappa)
names(w_GMRP_robust) <- colnames(X)
w_all_GMRP_robust_ellipsoid <- cbind(w_GMRP_robust)

dev.off()
# plot allocation
barplot(t(w_GMRP_robust), col = "red", legend = colnames(X), beside = TRUE,
        main = paste("Robust (ellipsoid) Global maximum return portfolio allocation, k=",kappa),
        xlab = "stocks", ylab = "dollars")


# Try multiple robust noisy solutions to check for sensitivity.
set.seed(357)
for (i in 1:6) {
  X_noisy <- rmvnorm(n = T, mean = mu, sigma = Sigma_SentInx)
  mu_noisy <- colMeans(X_noisy)
  Sigma_noisy <- cov(X_noisy)
  
  w_GMRP_robust_ellipsoid_noisy <- portfolioMaxReturnRobustEllipsoid(mu_noisy, Sigma_noisy, kappa)
  w_all_GMRP_robust_ellipsoid <- cbind(w_all_GMRP_robust_ellipsoid, w_GMRP_robust_ellipsoid_noisy)
}

# plot to compare the allocations
barplot(t(w_all_GMRP_robust_ellipsoid), col = 1:7, #legend = colnames(w_all_GMRP_robust_ellipsoid), 
        beside = TRUE, args.legend = list(bg = "white"),
        main = paste("Robust (ellipsoid) Global maximum return portfolio allocation, k=",kappa), xlab = "stocks", ylab = "dollars")


chart.StackedBar(t(w_all_GMRP_robust_ellipsoid), 
                 main = paste("Robust (ellipsoid) Global maximum return portfolio allocation, k=",kappa), 
                 ylab = "w", space = 0, border = NA)






