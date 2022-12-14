# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' @title sirscore
#' @description step 1 implemetation of algorithm sequentially and iteratively reweighted squares(SIRS).
#' @param A a matrix with dim n*p
#' @param y a vector with length n
#' @param a shape parameter of SICA penalty controling its maximum concavity 
#' @param x0 initial set
#' @param maxsize maximum size of the sparse model
#' @param delta param of SIRS
#' @param thresh threshold of the circulation. Default = 1e-6
#' @param maxiter maximum number of iterations. Default = 50
#' @param tol tolerance for updating.  Default = 1e-6
#' @return a coeffcient value in sparse recovery problem
#' @import Rcpp
#' @import RcppArmadillo
#' @useDynLib StatComp22010
#' @export
sirscore <- function(A, y, a, x0, maxsize, delta = 1e-6, thresh = 1e-6, maxiter = 50L, tol = 1e-6) {
    .Call('_StatComp22010_sirscore', PACKAGE = 'StatComp22010', A, y, a, x0, maxsize, delta, thresh, maxiter, tol)
}

#' @title sirs
#' @description Implemetation of algorithm sequentially and iteratively reweighted squares(SIRS).
#' @param A a matrix with dim n*p
#' @param y a vector with length n
#' @param x0 initial set
#' @param maxsize maximum size of the sparse model
#' @param eps scalar between 0 and 1, downweighting undesirable predictors
#' @param a shape parameter of SICA penalty controling its maximum concavity 
#' @param delta a small ridge parameter for stabilizing matrix inverse. Default = 1e-6
#' @param thresh threshold of the circulation. Default = 1e-6
#' @param maxiter maximum number of iterations. Default = 50
#' @param maxseq maximum number of sequential steps. Default = 50
#' @param tol tolerance for updating.  Default = 1e-6
#' @return a coeffcient value in sparse recovery problem
#' @examples
#' \dontrun{
#' library(MASS)
#' n = 100
#' p = 50
#' a = 0.4
#' rho = 0.5
#' mu = c(rep(0,p))
#' sigma1 <- matrix(0, p, p)
#' sigma1 <- rho ^ (abs(row(sigma1) - col(sigma1)))
#' A = mvrnorm(n,mu,sigma1)
#' beta = c(0.5,-0.5,1,-1.2,-1,rep(0,p-5))
#' y = A%*%beta
#' beta_hat=sirs(A,y,x0=rep(1,p),maxsize=min(ceiling(n/2),p),eps=1/p,a=a)
#' }
#' @import Rcpp
#' @import RcppArmadillo
#' @useDynLib StatComp22010
#' @export
sirs <- function(A, y, x0, maxsize, eps, a = 0.1, delta = 1e-6, thresh = 1e-6, maxiter = 50L, maxseq = 50L, tol = 1e-6) {
    .Call('_StatComp22010_sirs', PACKAGE = 'StatComp22010', A, y, x0, maxsize, eps, a, delta, thresh, maxiter, maxseq, tol)
}

#' @title Gibbs_binorm
#' @description Implemetation of Gibbs sampler on binorm situation.
#' @param N a matrix with dim n*p
#' @param sigma1 part of variance matrix Sigma, position [1,1]
#' @param sigma2 part of variance matrix Sigma, position [2,2]
#' @param mu1 part of mean vector, position [1]
#' @param mu2 part of mean vector, position [2]
#' @param rho part of variance matrix Sigma, position [1,2] or [2,1]
#' @param s1 s1 = sqrt(1-rho^2)*sigma1
#' @param s2 s2 = sqrt(1-rho^2)*sigma2
#' @param X initial matrix with first row given for the mean of giving start of Gibbs sampler  
#' @return a matrix of the result of Gibbs sampler of binorm situation
#' @examples
#' \dontrun{
#' #initialize constants and parameters
#' N = 5000 #length of chain
#' burn = 1000 #burn-in length
#' X = matrix(0, N, 2) #the chain, a bivariate sample
#' rho = 0.9 #correlation
#' mu1 = 0
#' mu2 = 0
#' sigma1 = 1
#' sigma2 = 1
#' s1 = sqrt(1-rho^2)*sigma1
#' s2 = sqrt(1-rho^2)*sigma2
#' # generate the chain
#' X[1, ] = c(mu1, mu2) #initialize
#' set.seed(1)
#' X = Gibbs_binorm(N,sigma1,sigma2,mu1,mu2,rho,s1,s2,X)
#' b = burn + 1
#' x = X[b:N, ]
#' plot(x, main="", cex=.5, xlab=bquote(X[t]),ylab=bquote(Y[t]), ylim=range(x[,2]))
#' }
#' @import Rcpp
#' @useDynLib StatComp22010
#' @export
Gibbs_binorm <- function(N, sigma1, sigma2, mu1, mu2, rho, s1, s2, X) {
    .Call('_StatComp22010_Gibbs_binorm', PACKAGE = 'StatComp22010', N, sigma1, sigma2, mu1, mu2, rho, s1, s2, X)
}

