## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(StatComp22010)

## -----------------------------------------------------------------------------
p=30
beta = c(c(1,-0.8,0.6,0,0,-1.5,-0.5,1.2),rep(0,p-8))
data1 = data_generate(n=50,p=30,beta = beta)

## -----------------------------------------------------------------------------
result = lrate_GIC(X=data1[[1]],y=data1[[2]])
learning_rate_GIC = result[1]
GIC_value = result[2]
#estimate beta
beta_est = estimate(X=data1[[1]],y=data1[[2]],learning_rate = unlist(learning_rate_GIC))
print(cbind(learning_rate_GIC,GIC_value))
print(cbind(beta,beta_est))

## -----------------------------------------------------------------------------
#compute estimation error.
prediction_error=pred_error(x=data1[[1]],y=data1[[2]],beta_hat = beta_est)
#compute estimation accuracy. 
estimation_accuracy=est_acc(beta_est,beta)
print(prediction_error)
print(estimation_accuracy)

## -----------------------------------------------------------------------------
library(MASS)
#generate the data for numeric study
n = 100
p = 50
a = 0.4
rho = 0.5
mu = c(rep(0,p))
sigma1 <- matrix(0, p, p)
sigma1 <- rho ^ (abs(row(sigma1) - col(sigma1)))
A = mvrnorm(n,mu,sigma1)
beta = c(0.5,-0.5,1,-1.2,-1,rep(0,p-5))
y = A%*%beta

#use function sirs
beta_hat=sirs(A,y,x0=rep(1,p),maxsize=min(ceiling(n/2),p),eps=1/p,a=a)
rbind(beta[1:5],beta_hat[1:5])
prediction_error = pred_error(exp(A),y,as.matrix(beta_hat))
estimation_accuracy = est_acc(beta_hat,beta)
print(prediction_error)
print(estimation_accuracy)

