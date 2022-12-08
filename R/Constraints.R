#' @title Algorithm to get a coefficient beta in compositional covariates problem
#' @name Algorithm1
#' @description Core pragram of the algorithm (Lin, W., Shi, P., Feng, R., Li, H, 2014), 
#' which solves the problem of variable selection in regression with compositional covariates.
#' @param x n x p design matrix with sample size n and dimensionality p
#' @param y n-vector of response
#' @param learning_rate regularization parameter of L1 penalty controling its maximum concavity
#' @param convergence_limit tolerance for updating.  Default = 1e-5
#' @return a coeffcient value in compositional covariates problem
#' @import wordspace
#' @useDynLib StatComp22010
#' @export
Algorithm1 <- function(x,y,learning_rate,convergence_limit = 1e-5){
  positive <- function(x){
    if (x>0) x
    else 0
  }
  
  soft_thresholding_operator <- function(t,learning_rate){
    sign(t)*positive((abs(t) - learning_rate))
  }
  
  step2_update <- function(n,p,j,z,y,beta,alpha,mu){
    su = 0
    for (i in 1:p){
      if (i != j) su = su + beta[i]*as.matrix(z[,i])
    }
    1/n*t(as.matrix(z[,j]))%*%(y-su)-mu*(sum(beta[-j]+alpha))
  }
  step2 <- function(n,p,z,y,mu,beta,alpha,convergence_limit,learning_rate){
    #judge if converge
    judge = numeric(p)
    beta0 = beta
    
    while (sum(judge)<p) {
      # for (j in 1:p) {
      #   beta_j_update = 1/((colNorms(z)[j])^2/n+mu)*soft_thresholding_operator(step2_update(n,p,j,z,y,beta0,alpha,mu),learning_rate)
      #   if (abs(beta_j_update-beta0[j])>convergence_limit) beta[j] = beta_j_update
      #   else judge[j] = 1
      # beta0 = beta
      # }
      for (j in 1:p) {
        beta_j_update = 1/((wordspace::colNorms(z)[j])^2/n+mu)*soft_thresholding_operator(step2_update(n,p,j,z,y,beta,alpha,mu),learning_rate)
        if (abs(beta_j_update-beta[j])>convergence_limit) beta[j] = beta_j_update
        else judge[j] = 1
      }
    }
    return(beta)
  }
  #need matrix as input X,y
  
  #Initialize
  n = dim(x)[1]
  p = dim(x)[2]
  z = log(x)
  judge = 1
  beta = numeric(p);alpha = 0;mu = 1;k = 0;j = 1
  
  while (judge) {
    
    
    beta_update = step2(n,p,z,y,mu,beta,alpha,convergence_limit,learning_rate)
    if (sum((beta_update-beta)^2)>convergence_limit){
      #alpha update
      beta = beta_update
      alpha = alpha+sum(beta)
      k = k+1 
      #print(k)
      #print(alpha)
      #print(beta)
    }
    else {
      judge = 0
      beta = beta_update  
    }
  }
  return(beta)
}

#' @title GIC computation.
#' @name GIC
#' @description GIC value computation, which is defined by (Lin, W., Shi, P., Feng, R., Li, H, 2014).
#' @param lambda regularization parameter of L1 penalty
#' @param z n x p log design matrix with sample size n and dimensionality p
#' @param y n-vector of response
#' @param beta_lambda  coeffcient value of regression under regularization parameter lambda
#' @return a GIC value under the data given
#' @examples 
#' \dontrun{
#' p=30
#' beta = c(c(1,-0.8,0.6,0,0,-1.5,-0.5,1.2),rep(0,p-8))
#' data1 = data_generate(n=50,p=30,beta = beta)
#' result = lrate_GIC(X=data1[[1]],y=data1[[2]])
#' learning_rate_GIC = result[1]
#' GIC_value = result[2]
#' 
#' #estimate beta
#' beta_est = estimate(X=data1[[1]],y=data1[[2]],learning_rate = unlist(learning_rate_GIC))
#' 
#' #compute estimation error.
#' pred_error(x=data1[[1]],y=data1[[2]],beta_hat = beta_est)
#' 
#' #compute estimation accuracy. 
#' est_acc(beta_est,beta)
#' }
#' @import wordspace
#' @useDynLib StatComp22010
#' @export
GIC <- function(lambda,z,y,beta_lambda){
  #select lambda by min GIC
  n = dim(z)[1]
  p = dim(z)[2]
  #need matrix as input:y,z
  wordspace::colNorms(y-z%*%beta_lambda)^2/n+(sum(beta_lambda!=0)-1)*log(log(n))/n*log(max(p,n))
}

#' @title Select regularization parameter under minimum GIC
#' @name lrate_GIC
#' @description by minimum GIC value to select the best regularization parameter.
#' @param X n x p design matrix with sample size n and dimensionality p
#' @param y n-vector of response
#' @param tol  tolerance for updating.  Default = 1e-6 
#' @return a list including the beat regularization parameter selected and the minimum GIC reached
#' @examples 
#' \dontrun{
#' p=30
#' beta = c(c(1,-0.8,0.6,0,0,-1.5,-0.5,1.2),rep(0,p-8))
#' data1 = data_generate(n=50,p=30,beta = beta)
#' result = lrate_GIC(X=data1[[1]],y=data1[[2]])
#' learning_rate_GIC = result[1]
#' GIC_value = result[2]
#' 
#' #estimate beta
#' beta_est = estimate(X=data1[[1]],y=data1[[2]],learning_rate = unlist(learning_rate_GIC))
#' 
#' #compute estimation error.
#' pred_error(x=data1[[1]],y=data1[[2]],beta_hat = beta_est)
#' 
#' #compute estimation accuracy. 
#' est_acc(beta_est,beta)
#' }
#' @import stats
#' @useDynLib StatComp22010
#' @export
lrate_GIC <- function(X,y,tol=1e-6){
  result = stats::optimize(function(x) GIC(x,log(X),y,Algorithm1(x=X,y=y,learning_rate=x,convergence_limit = tol)),lower = 0,upper = 1,maximum = F) 
  return(list(learning_rate_GIC=result$minimum,GIC=result$objective))
}

#' @title Algorithm to get a coefficient beta in compositional covariates problem
#' @name estimate
#' @description  Solves the problem of variable selection in regression with compositional covariates, and estimate the cofficient beta.
#' @param X n x p design matrix with sample size n and dimensionality p
#' @param y n-vector of response
#' @param learning_rate regularization parameter of L1 penalty controling its maximum concavity
#' @param tol tolerance for updating.  Default = 1e-5
#' @return a coefficient value in compositional covariates problem
#' @examples 
#' \dontrun{
#' p=30
#' beta = c(c(1,-0.8,0.6,0,0,-1.5,-0.5,1.2),rep(0,p-8))
#' data1 = data_generate(n=50,p=30,beta = beta)
#' result = lrate_GIC(X=data1[[1]],y=data1[[2]])
#' learning_rate_GIC = result[1]
#' GIC_value = result[2]
#' 
#' #estimate beta
#' beta_est = estimate(X=data1[[1]],y=data1[[2]],learning_rate = unlist(learning_rate_GIC))
#' 
#' #compute estimation error.
#' pred_error(x=data1[[1]],y=data1[[2]],beta_hat = beta_est)
#' 
#' #compute estimation accuracy. 
#' est_acc(beta_est,beta)
#' }
#' @useDynLib StatComp22010
#' @export
estimate <- function(X,y,learning_rate,tol = 1e-5) {
  Algorithm1(x=X,y=y,learning_rate = learning_rate,convergence_limit = tol)
}

#' @title Generate the data for numeric study
#' @name data_generate
#' @description  generate the data for numeric study in the exact paper and for example.
#' @param n sample size n
#' @param p dimension
#' @param beta true coefficient beta
#' @param rho for generate covariance matrix. Default =0.2
#' @param sig covariance for white noise. Default = 0.5
#' @return a list of data generated including data X and response y in the form of matrix
#' @examples 
#' \dontrun{
#' p=30
#' beta = c(c(1,-0.8,0.6,0,0,-1.5,-0.5,1.2),rep(0,p-8))
#' data1 = data_generate(n=50,p=30,beta = beta)
#' result = lrate_GIC(X=data1[[1]],y=data1[[2]])
#' learning_rate_GIC = result[1]
#' GIC_value = result[2]
#' 
#' #estimate beta
#' beta_est = estimate(X=data1[[1]],y=data1[[2]],learning_rate = unlist(learning_rate_GIC))
#' 
#' #compute estimation error.
#' pred_error(x=data1[[1]],y=data1[[2]],beta_hat = beta_est)
#' 
#' #compute estimation accuracy. 
#' est_acc(beta_est,beta)
#' }
#' @import MASS
#' @import stats
#' @useDynLib StatComp22010
#' @export
data_generate <- function(n,p,beta,rho=0.2,sig=0.5){
  #judge p>5?
  if (p<=5) theta = rep(log(0.5*p),p)
  else theta = c(rep(log(0.5*p),5),rep(0,p-5))
  
  #initial
  sigma = matrix(0,p,p)
  
  for (i in 1:p){
    for (j in 1:p) {
      sigma[i,j] = rho^(abs(i-j))
    }
  }
  
  W = MASS::mvrnorm(n,mu = theta,Sigma = sigma)
  X = apply(exp(W), 2, function(x) x/sum(x))
  #X = t(apply(exp(W), 1, function(x) x/sum(x)))
  #need p>8
  #beta_star = c(c(1,-0.8,0.6,0,0,-1.5,-0.5,1.2),rep(0,p-8))
  Z = log(X)
  y = Z%*%beta+stats::rnorm(n,0,sig)
  return(list(X,y))
}

#' @title Evaluate the performance of estimation
#' @name pred_error
#' @description  compute the prediction error of estimation.
#' @param x n x p design matrix with sample size n and dimensionality p
#' @param y n-vector of response
#' @param beta_hat estimated coefficient beta
#' @return the predict error
#' @examples 
#' \dontrun{
#' p=30
#' beta = c(c(1,-0.8,0.6,0,0,-1.5,-0.5,1.2),rep(0,p-8))
#' data1 = data_generate(n=50,p=30,beta = beta)
#' result = lrate_GIC(X=data1[[1]],y=data1[[2]])
#' learning_rate_GIC = result[1]
#' GIC_value = result[2]
#' 
#' #estimate beta
#' beta_est = estimate(X=data1[[1]],y=data1[[2]],learning_rate = unlist(learning_rate_GIC))
#' 
#' #compute estimation error.
#' pred_error(x=data1[[1]],y=data1[[2]],beta_hat = beta_est)
#' 
#' #compute estimation accuracy. 
#' est_acc(beta_est,beta)
#' }
#' @import wordspace
#' @useDynLib StatComp22010
#' @export
pred_error <- function(x,y,beta_hat){
  #need matrix input :y,x,beta
  wordspace::colNorms(y-log(x)%*%beta_hat)^2/dim(x)[1]
}

#' @title Evaluate the performance of estimation
#' @name est_acc
#' @description  compute the accuracy of estimation.
#' @param beta_hat estimated coefficient beta
#' @param beta_true true coefficient beta
#' @return a list of loss_1,loss_2 and loss_infinity under different norms
#' @examples 
#' \dontrun{
#' p=30
#' beta = c(c(1,-0.8,0.6,0,0,-1.5,-0.5,1.2),rep(0,p-8))
#' data1 = data_generate(n=50,p=30,beta = beta)
#' result = lrate_GIC(X=data1[[1]],y=data1[[2]])
#' learning_rate_GIC = result[1]
#' GIC_value = result[2]
#' 
#' #estimate beta
#' beta_est = estimate(X=data1[[1]],y=data1[[2]],learning_rate = unlist(learning_rate_GIC))
#' 
#' #compute estimation error.
#' pred_error(x=data1[[1]],y=data1[[2]],beta_hat = beta_est)
#' 
#' #compute estimation accuracy. 
#' est_acc(beta_est,beta)
#' }
#' @import wordspace
#' @useDynLib StatComp22010
#' @export
est_acc <- function(beta_hat,beta_true){
  #l_q losses:q=1,2,inifty
  beta_star = beta_true
  loss_1 = wordspace::colNorms(as.matrix(beta_hat-beta_star),method = "manhattan",p=1)
  loss_2 = wordspace::colNorms(as.matrix(beta_hat-beta_star),method = "euclidean",p=2)
  loss_infty = wordspace::colNorms(as.matrix(beta_hat-beta_star),method = "maximum",p=Inf)
  return(list(loss_1=loss_1,loss2=loss_2,loss_infty=loss_infty))
}