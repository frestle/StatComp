% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Constraints.R
\name{estimate}
\alias{estimate}
\title{Algorithm to get a coefficient beta in compositional covariates problem}
\usage{
estimate(X, y, learning_rate, tol = 1e-05)
}
\arguments{
\item{X}{n x p design matrix with sample size n and dimensionality p}

\item{y}{n-vector of response}

\item{learning_rate}{regularization parameter of L1 penalty controling its maximum concavity}

\item{tol}{tolerance for updating.  Default = 1e-5}
}
\value{
a coefficient value in compositional covariates problem
}
\description{
Solves the problem of variable selection in regression with compositional covariates, and estimate the cofficient beta.
}
\examples{
\dontrun{
p=30
beta = c(c(1,-0.8,0.6,0,0,-1.5,-0.5,1.2),rep(0,p-8))
data1 = data_generate(n=50,p=30,beta = beta)
result = lrate_GIC(X=data1[[1]],y=data1[[2]])
learning_rate_GIC = result[1]
GIC_value = result[2]

#estimate beta
beta_est = estimate(X=data1[[1]],y=data1[[2]],learning_rate = unlist(learning_rate_GIC))

#compute estimation error.
pred_error(x=data1[[1]],y=data1[[2]],beta_hat = beta_est)

#compute estimation accuracy. 
est_acc(beta_est,beta)
}
}
