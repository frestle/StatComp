\name{StatComp22010-package}
\alias{StatComp22010-package}
\alias{StatComp22010}
\docType{package}
\title{
This package provides implementation of one method of variable selection in regression with compositional covariates, sequentially and iteratively reweighted squares(SIRS) with Rcpp and all of the homework of the cource "statistic computing".
}
\description{

 With a method raised by WEI LIN, PIXU SHI, RUI FENG AND HONGZHE LI, this package provide a function estimate() to do regression with compositional covariates so as to estimate the coefficients, note that before using this function, matrix X and vector y should be centralized first; moreover, this package provides a function sirs() written by Rcpp to implement SIRS algorithm, in addition, all of my homework of the cource "statistic computing" is included. 
}
\details{ 

estimate() requires input of matrix X and vector y and learning rate, and will estimate the coffcient beta that satisfies X*beta=y under the learning rate, function lrate_GIC() will search the best learning rate that minimize the GIC value, requires the input of matrix X and vector y, and return the minimum GIC and the corresponding learning rate. function sirs() implement the algorithm SIRS, require the input matrix X, vector y and parameter a for the algorithm, return the estimate of beta similar to function estimate(),in addition, there're some functions to evaluate the performance of estimation.  

}
\author{
22010 <zhangxiao@mail.ustc.edu.cn>
}
\examples{
 \dontrun{

#function1
# initial dims
p=30
beta = c(c(1,-0.8,0.6,0,0,-1.5,-0.5,1.2),rep(0,p-8))

#generate data

data1 = data_generate(n=50,p=30,beta = beta)
result = lrate_GIC(X=data1[[1]],y=data1[[2]])
learning_rate_GIC = result[1]
GIC_value = result[2]

# estimate beta

beta_est = estimate(X=data1[[1]],y=data1[[2]],learning_rate = unlist(learning_rate_GIC))

#compute estimation error.

pred_error(x=data1[[1]],y=data1[[2]],beta_hat = beta_est)

# compute estimation accuracy. 

est_acc(beta_est,beta)


#function2

library(MASS)

#initial dims

n = 100
p = 50
a = 0.4
rho = 0.5
mu = c(rep(0,p))
sigma1 <- matrix(0, p, p)
sigma1 <- rho ^ (abs(row(sigma1) - col(sigma1)))

# generate data

A = mvrnorm(n,mu,sigma1)
beta = c(0.5,-0.5,1,-1.2,-1,rep(0,p-5))
y = A %*% beta

beta_hat=sirs(A,y,x0=rep(1,p),maxsize=min(ceiling(n/2),p),eps=1/p,a=a)
}
}

\references{
\itemize{

 \item [1] Lv J, Fan Y. A unified approach to model selection and sparse recovery using regularized least squares. Ann. Statist. 2009;37:3498–3528.

 \item [2] Lin, W., Shi, P., Feng, R., Li, H. Variable selection in regression with compositional covariates. Biometrika. 2014;101:785–797.
}
}