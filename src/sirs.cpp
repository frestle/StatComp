#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

//' @title sirscore
//' @description step 1 implemetation of algorithm sequentially and iteratively reweighted squares(SIRS).
//' @param A a matrix with dim n*p
//' @param y a vector with length n
//' @param a shape parameter of SICA penalty controling its maximum concavity 
//' @param x0 initial set
//' @param maxsize maximum size of the sparse model
//' @param delta param of SIRS
//' @param thresh threshold of the circulation. Default = 1e-6
//' @param maxiter maximum number of iterations. Default = 50
//' @param tol tolerance for updating.  Default = 1e-6
//' @return a coeffcient value in sparse recovery problem
//' @import Rcpp
//' @import RcppArmadillo
//' @useDynLib StatComp22010
//' @export
// [[Rcpp::export]]
arma::vec sirscore(arma::mat A,arma::vec y,double a,arma::vec x0,int maxsize,double delta=1e-6,double thresh=1e-6, int maxiter=50, double tol=1e-6) {
  int n = A.n_rows;
  int p = A.n_cols;
  if (x0.n_elem==0){
    x0=arma::ones(p);
  }
  
  arma::vec x = x0;
  arma::mat D = arma::eye(p,p);
  arma::mat D1=D;
  for (int i=0;i<p;i++){
    D(i,i)=std::abs(x(i))*(a + std::abs(x(i)))/(a + 1);
  }
  int k = 1;
  int update = 1;
  while (update > tol && k <= maxiter){
    k = k + 1;
    arma::vec xold = x;
    double nlog = std::log(n);
    if (p <= std::ceil(nlog)*n){
      D1 = arma::sqrt(D);
      x = D1 * arma::inv(delta*arma::eye(p,p) + D1 * A.t() * A * D1) * D1 * A.t() * y;
    }else {
      x = D * A.t() * arma::inv(delta*arma::eye(n,n) + A * D * A.t()) * y; 
    }
    
    update = arma::sum(arma::pow(x - xold,2));
    for (int i=0;i<p;i++){
      double xi=x(i);
      D(i,i)=std::abs(xi)*(a + std::abs(xi))/(a + 1);
    }
  }
  
  arma::uvec thre = arma::find(arma::abs(x) > thresh);
  
  int pp = x.n_elem;
  for (int i=0;i<pp;i++){
    if (arma::any(thre-i) == false){
      x(i) = 0;
    }
  }
  int num = thre.n_elem;
  arma::vec xp;
  
  if (num <= maxsize){
    xp = x;
    arma::uvec estmod = arma::find(xp!=0);
    int len = estmod.n_elem;
    arma::mat A_mod = A.cols(estmod);
    arma::vec xpp = arma::inv(A_mod.t() * A_mod) * A_mod.t() * y;
    for (int i=0;i<len;i++){
      xp(estmod(i)) = xpp(i);
    }
    arma::uvec xpthre = arma::find(arma::abs(xp) > thresh);
    for (int i=0;i<pp;i++){
      if (arma::any(thre-i) == false){
        x(i) = 0;
      }
    }
  }else{ 
    xp = x;
  }
  return xp;
}


//' @title sirs
//' @description Implemetation of algorithm sequentially and iteratively reweighted squares(SIRS).
//' @param A a matrix with dim n*p
//' @param y a vector with length n
//' @param x0 initial set
//' @param maxsize maximum size of the sparse model
//' @param eps scalar between 0 and 1, downweighting undesirable predictors
//' @param a shape parameter of SICA penalty controling its maximum concavity 
//' @param delta a small ridge parameter for stabilizing matrix inverse. Default = 1e-6
//' @param thresh threshold of the circulation. Default = 1e-6
//' @param maxiter maximum number of iterations. Default = 50
//' @param maxseq maximum number of sequential steps. Default = 50
//' @param tol tolerance for updating.  Default = 1e-6
//' @return a coeffcient value in sparse recovery problem
//' @examples
//' \dontrun{
//' library(MASS)
//' n = 100
//' p = 50
//' a = 0.4
//' rho = 0.5
//' mu = c(rep(0,p))
//' sigma1 <- matrix(0, p, p)
//' sigma1 <- rho ^ (abs(row(sigma1) - col(sigma1)))
//' A = mvrnorm(n,mu,sigma1)
//' beta = c(0.5,-0.5,1,-1.2,-1,rep(0,p-5))
//' y = A%*%beta
//' beta_hat=sirs(A,y,x0=rep(1,p),maxsize=min(ceiling(n/2),p),eps=1/p,a=a)
//' }
//' @import Rcpp
//' @import RcppArmadillo
//' @useDynLib StatComp22010
//' @export
// [[Rcpp::export]]
arma::vec sirs(arma::mat A,arma::vec y,arma::vec x0,int maxsize,double eps,double a=0.1,double delta=1e-6,double thresh=1e-6,int maxiter=50,int maxseq=50,double tol=1e-6){
  //int n=A.n_rows;
  int p=A.n_cols;
  
  //vec ma = linspace(0, 1, 2);
  //ma(0) = ceil(n/2);
  //ma(1)=p;
  //maxsize=min(ma);
  
  //eps=1/p;
  //x0=ones(p);
  arma::vec x_ini = x0;
  int rep = 1;
  int move = 1;
  arma::vec xp0 = arma::zeros(p);
  arma::vec xp;
  while (move > tol && rep <= maxseq){
    xp = sirscore(A, y, a, x_ini, maxsize, delta, thresh, maxiter, tol);
    arma::uvec roo = arma::find(xp!=0);
    int num = roo.n_elem;
    if (num <= maxsize){
      break;
    }else{
      arma::uvec estmod = arma::find(xp!=0);
      int len = estmod.n_elem;
      arma::vec xd=arma::zeros(len);
      for(int i=0;i<len;i++){
        xd(i)=xp(estmod(i));
      }
      arma::vec xd_abs = arma::abs(xd);
      arma::vec kth_vec = arma::sort(xd_abs);
      double kth = kth_vec(rep-1);
      
      x_ini = eps*arma::ones(p);
      for(int i=0;i<len;i++){
        if(xd(i)>=kth){
          x_ini(i)=1;
        }
      }
      move = arma::sum(arma::pow(xp - xp0,2));
      xp0 = xp;
      rep = rep + 1;
      
    }
  }
  return xp;
}

//' @title Gibbs_binorm
//' @description Implemetation of Gibbs sampler on binorm situation.
//' @param N a matrix with dim n*p
//' @param sigma1 part of variance matrix Sigma, position [1,1]
//' @param sigma2 part of variance matrix Sigma, position [2,2]
//' @param mu1 part of mean vector, position [1]
//' @param mu2 part of mean vector, position [2]
//' @param rho part of variance matrix Sigma, position [1,2] or [2,1]
//' @param s1 s1 = sqrt(1-rho^2)*sigma1
//' @param s2 s2 = sqrt(1-rho^2)*sigma2
//' @param X initial matrix with first row given for the mean of giving start of Gibbs sampler  
//' @return a matrix of the result of Gibbs sampler of binorm situation
//' @examples
//' \dontrun{
//' #initialize constants and parameters
//' N = 5000 #length of chain
//' burn = 1000 #burn-in length
//' X = matrix(0, N, 2) #the chain, a bivariate sample
//' rho = 0.9 #correlation
//' mu1 = 0
//' mu2 = 0
//' sigma1 = 1
//' sigma2 = 1
//' s1 = sqrt(1-rho^2)*sigma1
//' s2 = sqrt(1-rho^2)*sigma2
//' # generate the chain
//' X[1, ] = c(mu1, mu2) #initialize
//' set.seed(1)
//' X = Gibbs_binorm(N,sigma1,sigma2,mu1,mu2,rho,s1,s2,X)
//' b = burn + 1
//' x = X[b:N, ]
//' plot(x, main="", cex=.5, xlab=bquote(X[t]),ylab=bquote(Y[t]), ylim=range(x[,2]))
//' }
//' @import Rcpp
//' @useDynLib StatComp22010
//' @export
// [[Rcpp::export]]
NumericMatrix Gibbs_binorm(int N, double sigma1, double sigma2, double mu1, double mu2, double rho, double s1, double s2, NumericMatrix X) {
  
  for (int i=1;i<N;i++) {
    double x2 = X(i-2, 1);
    double m1 = mu1 + rho * (x2 - mu2) * sigma1/sigma2;
    X(i-1, 0) = rnorm(1, m1, s1)[0];
    double x1 = X(i-1, 0);
    double m2 = mu2 + rho * (x1 - mu1) * sigma2/sigma1;
    X(i-1, 1) = rnorm(1, m2, s2)[0];
  }
  return X;
}
