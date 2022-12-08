## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(StatComp22010)

## -----------------------------------------------------------------------------
options (warn = -1)
library("knitr")
dataset1 <- mtcars

## -----------------------------------------------------------------------------
row.names(dataset1)

## -----------------------------------------------------------------------------
disp <- mtcars$disp 
hist_disp <- hist(disp, col="blue", breaks=9, xlab="disp", main="histgram of disp in dataset\"mtcars\"")

## -----------------------------------------------------------------------------
kable(head(dataset1))

## -----------------------------------------------------------------------------
n = 1e3 #n for all problems thus without distinction.
a_pareto_Q1 = 2
b_pareto_Q1 = 2
set.seed(91501)
u_Q1 = runif(n)
x_Q1 = b_pareto_Q1*(1-u_Q1)^(-1/a_pareto_Q1) #we know it from step1.

## -----------------------------------------------------------------------------
hist(x_Q1, prob=TRUE, main="Inverse Transform Method for Pareto(2,2)",xlab="x") #we know the function from step1.
y_Q1 = seq(b_pareto_Q1, 100, .01)
lines(y_Q1, a_pareto_Q1*b_pareto_Q1^a_pareto_Q1*(1/y_Q1)^(a_pareto_Q1+1),lwd=1.5) 

## -----------------------------------------------------------------------------
Gen_beta = function(n,a,b){
  k = 0
  y = numeric(n)
  while (k < n){
    u = runif(1)
    x = runif(1) #the envelope distribution.
    #use optimize to find the best c in acceptance-rejection method with accuracy 0.001.
    max_beta = optimize(f<-function(x) {x^(a-1)*(1-x)^(b-1)},interval=c(0,1),tol=0.001,maximum = TRUE)$maximum 
    if (x^(a-1)*(1-x)^(b-1)>(max_beta+0.001)*u){
      #accept x
      k = k + 1
      y[k] = x
    }
  }
  return(y)
}

## -----------------------------------------------------------------------------
a_Q2 = 3
b_Q2 = 2
set.seed(91802)
sample_Q2 = Gen_beta(n,a_Q2,b_Q2)

## -----------------------------------------------------------------------------
hist(sample_Q2, prob=TRUE, main="Acceptance-Rejection Method for Beta(3,2)",xlab="x") #we know the function from step3.
y_Q2 = seq(0, 1, .01)
lines(y_Q2, 12*y_Q2^2*(1-y_Q2),lwd=1.5) 

## -----------------------------------------------------------------------------
r_Q3 = 4
beta_Q3 = 2
set.seed(91503)
lambda_Q3 = rgamma(n,shape=4,rate=beta_Q3)
y_Q3 = rexp(n,lambda_Q3)

## -----------------------------------------------------------------------------
set.seed(91504)
lambda_Q4 = rgamma(n,shape=r_Q3,rate=beta_Q3)
y_Q4 = rexp(n,lambda_Q4)
hist(y_Q4, prob=TRUE, main="Empirical and Theoretical Pareto Distributions",xlab="x")
z_Q4 = seq(0, 100, .01)
lines(z_Q4, r_Q3*beta_Q3^r_Q3*(1/(z_Q4+beta_Q3))^(r_Q3+1),lwd=1.5) #we know the function from step3.

## -----------------------------------------------------------------------------
options (warn = -1)
library(knitr)
m = c(1e4,2*1e4,4*1e4,6*1e4,8*1e4) #numbers set
n = 100 #times set
quick_sort<-function(x){
  num = length(x)
  if(num == 0||num == 1){return(x)
  }else{
    a = x[1]
    y = x[-1]
    lower = y[y<a]
    upper = y[y>=a]
    return(c(quick_sort(lower),a,quick_sort(upper)))} #form a loop statement
}
quick_sort_1<-function(x){ #prepare for the calculation of computation time
  y = sample(1:x) #for the step2 operation seed is not set
  return(system.time(quick_sort(y))[1])
}
time_Q1 = lapply(m,quick_sort_1)
result_Q1 = data.frame(matrix(c(m,time_Q1),ncol = 2)) #show the result in the form of table
colnames(result_Q1) = c("number","time(s)")
kable(result_Q1,caption = "Fast Sorting Algorithm Time for Randomly Permuted Numbers")

## -----------------------------------------------------------------------------
m_loop = rep(m,n)
time_Q1_loop = unlist(lapply(m_loop,quick_sort_1)) #just for habit,avoid using loops with "for","while",etc.
time_Q1_aver = rowMeans(matrix(time_Q1_loop,nrow = 5))
result_Q1_aver = data.frame(matrix(c(m,time_Q1_aver),ncol = 2)) #show the result in the form of table
colnames(result_Q1_aver) = c("number","aver_time over 100(s)")
kable(result_Q1_aver,caption = "Computation Time Averaged over 100 Simulations")

## -----------------------------------------------------------------------------
num_trans = m*log(m)
library(ggplot2)
data_Q1 = data.frame(num_trans,time_Q1_aver)
p <- ggplot(data_Q1,aes(x=num_trans,y=time_Q1_aver)) + geom_point(shape=19) + xlab("nlog(n)") + ylab("time")+geom_smooth(method = lm)
p+ggtitle("Regress time on nlog(n)")+theme(plot.title = element_text(hjust = 0.5))

## -----------------------------------------------------------------------------
e = exp(1)
100*(e^2-3*e+1)/(0.5*(e^2-1)-(e-1)^2)

## -----------------------------------------------------------------------------
n_Q3 = 1e4
set.seed(9233)
MC.Phi <- function(R = 1e4, antithetic = TRUE){ #two Ways to estimate theta
  u = runif(R/2)
  if (!antithetic) v = runif(R/2) else v = 1 - u
  u = c(u, v)
  g = exp(u)
  mean(g)
}
MC1 = MC2 = numeric(n_Q3)
for (i in 1:n_Q3) {
  MC1[i] = MC.Phi(R = n_Q3, anti = FALSE)
  MC2[i] = MC.Phi(R = n_Q3)
}
result_Q3 = data.frame(matrix(c(mean(MC1),mean(MC2)),ncol = 2)) #show the result in the form of table
colnames(result_Q3) = c("Simple Monte Carlo","Antithetic variate")
rownames(result_Q3) = "theta_hat"
kable(result_Q3,caption = "Value of Theta Estimated in Different Ways")

## -----------------------------------------------------------------------------
percent = 100*(var(MC1)-var(MC2))/var(MC1)
percent

## -----------------------------------------------------------------------------
options(warn = -1)
library(ggplot2)
library(reshape2)
library(VGAM)
library(knitr)
#define target function g(x).
G_Q1 <- function(x){ 
  x^2/sqrt(2*pi)*exp(-x^2/2)*(x>1)
}
x = seq(1,10,0.01)
#we choose importance functions:chisq(4) and rayleigh(sqrt(2))
data_Q1 = data.frame(x,g=G_Q1(x),chisq=dchisq(x,df=4),
                     rayleigh=drayleigh(x,scale = sqrt(2),log=F))
#codes below are to graph the funtions.
data.plot = melt(data_Q1,id="x")
colnames(data.plot) <- c("x","func","value")
ggplot(data = data.plot,aes(x=x,y=value,group = func,color=func))+
  geom_line(lwd=2)+
  xlab("x")+
  ylab("value")+
  theme_bw()+
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        legend.position = c(.775,.815),
        legend.box.background = element_rect(color="black"))+
  scale_x_continuous(limits = c(1,10),breaks = seq(1,10,1))+
  ggtitle("Graph of Imporance Functions and Target Function")+
  theme(plot.title = element_text(hjust = 0.5))

## -----------------------------------------------------------------------------
m = 1e4 
set.seed(0930)
x1 = rchisq(m,df=4)
x2 = rrayleigh(m,scale=sqrt(2))
mean_Q1 = round(c(mean(G_Q1(x1)/dchisq(x1,df=4)),
                  mean(G_Q1(x2)/drayleigh(x2,scale=sqrt(2)))),3)
variance_Q1 = round(c(var(G_Q1(x1)/dchisq(x1,df=4)),
                      var(G_Q1(x2)/drayleigh(x2,scale=sqrt(2)))),3)
#codes below are to show the result in the form of table with knitr.
result_Q1 = data.frame(mean_Q1,variance_Q1)
colnames(result_Q1) = c("mean","var")
rownames(result_Q1) = c("chisq","rayleigh")
kable(result_Q1,caption = "Result of Importance Sampling")

## -----------------------------------------------------------------------------
#change the data and plot again like step1.
data_Q1.2 = data.frame(x,chisq=G_Q1(x)/dchisq(x,df=4),
                       rayleigh=G_Q1(x)/drayleigh(x, scale = sqrt(2),log=F))
data.plot.2 = melt(data_Q1.2,id="x")
colnames(data.plot.2) <- c("x","func","value")
ggplot(data = data.plot.2,aes(x=x,y=value,group = func,color=func))+
  geom_line(lwd=2)+
  xlab("x")+
  ylab("value")+
  theme_bw() +
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        legend.position = c(.775,.815),
        legend.box.background = element_rect(color="black"))+
  scale_x_continuous(limits = c(1,10),breaks = seq(1,10,1))+
  ggtitle("Ratio Function of Imporance Functions and Target Function")+
  theme(plot.title = element_text(hjust = 0.5))

## -----------------------------------------------------------------------------
#codes below are copied from Statistical Computing with R, page141.
set.seed(0930)
theta.hat <- se <- numeric(5)
g <- function(x) {
exp(-x - log(1+x^2)) * (x > 0) * (x < 1)
}
x <- runif(m) #using f0
fg <- g(x)
theta.hat[1] <- mean(fg)
se[1] <- sd(fg)
x <- rexp(m, 1) #using f1
fg <- g(x) / exp(-x)
theta.hat[2] <- mean(fg)
se[2] <- sd(fg)
x <- rcauchy(m) #using f2
i <- c(which(x > 1), which(x < 0))
x[i] <- 2 #to catch overflow errors in g(x)
fg <- g(x) / dcauchy(x)
theta.hat[3] <- mean(fg)
se[3] <- sd(fg)
u <- runif(m) #f3, inverse transform method
x <- - log(1 - u * (1 - exp(-1)))
fg <- g(x) / (exp(-x) / (1 - exp(-1)))
theta.hat[4] <- mean(fg)
se[4] <- sd(fg)
u <- runif(m) #f4, inverse transform method
x <- tan(pi * u / 4)
fg <- g(x) / (4 / ((1 + x^2) * pi))
theta.hat[5] <- mean(fg)
se[5] <- sd(fg)
result_Q2.1 = round(rbind(theta.hat, se),3)
colnames(result_Q2.1) = c("f_0","f_1","f_2","f_3","f_4")
kable(result_Q2.1,
      caption = "Result of Importance Sampling")

## -----------------------------------------------------------------------------
set.seed(0930)
f.Q2 <- function(x){ #Density f_j(x) = k*f(x) where k=5.
  5*exp(-x)/(1-exp(-1))*(x<1)*(x>0)
}
 F.revers1 <- function(x,a){ #Define the reverse pmf of the density above.
   -log(exp(-a)-0.2*x*(1-exp(-1)))
 }
 F.revers <- function(x){ #Define the reverse pmf of f_3 to split the interval.
   -log(1-x*(1-exp(-1)))
 }
inteval = unlist(c(0,lapply(1:4/5, F.revers),1))
inteval
theta.hat.str=va=NULL
for (i in 1:5) {
  u = runif(m/5)            
  x = F.revers1(u,inteval[i]) #Generate random numbers by f_j.
  fg = g(x)*(x>=inteval[i])*(x<inteval[i+1])/f.Q2(x) 
  theta.hat.str[i] = mean(fg)
  va[i] = var(fg)
}
result_Q2.2 = round(data.frame(sum(theta.hat.str),sqrt(sum(va))),4)
colnames(result_Q2.2) = c("Mean","Se")
kable(result_Q2.2,
      caption = "Result of Stratified Importance Sampling")

## -----------------------------------------------------------------------------
result_Q2 = matrix(c(sum(theta.hat.str),theta.hat[4],sqrt(sum(va)),se[4]),
                   ncol = 2,dimnames=list(c("EX5.13","EX5.10"),c("Mean","Se")))
kable(round(result_Q2,4),
      caption = "Comparison between the Results of 2 Examples")

## -----------------------------------------------------------------------------
#options(warn = -1)
library(knitr)
set.seed(109)
n = 1000
alpha = .05
UCL <- replicate(1000, expr = {
  x = rlnorm(n,0,1) #set parameter meanlog = 0, sdlog = 1.
   sqrt(n)*abs(mean(log(x)))/var(log(x))/qt(alpha/2, df = n-1, lower.tail = F)
} )
mean(UCL<=1)

## -----------------------------------------------------------------------------
count5test <- function(x, y) {
X <- x - mean(x)
Y <- y - mean(y)
outx <- sum(X > max(Y)) + sum(X < min(Y))
outy <- sum(Y > max(X)) + sum(Y < min(X))
# return 1 (reject) or 0 (do not reject H0)
return(as.integer(max(c(outx, outy)) > 5))
}
#define the function to give the result of different tests and sample sizes.
Q2test <- function(n){
 # generate samples under H1 to estimate power
 m = 1000
 sigma1 = 1
 sigma2 = 1.5
 alpha = 0.055
 power = replicate(m, expr={
 x = rnorm(n, 0, sigma1)
 y = rnorm(n, 0, sigma2)
 c(count5test(x,y),as.integer(var.test(x,y)$p.value <= alpha))
 })
 return(c(mean(power[1,]),mean(power[2,])))
}

## -----------------------------------------------------------------------------
set.seed(109)
resultQ2 = lapply(c(10,100,1000),Q2test)
#to show the result in the form of table.
re_t = matrix(unlist(resultQ2),nrow = 2)
rownames(re_t) = c("Count Five Test","F test")
colnames(re_t) = c(10,100,1000)
kable(re_t, caption = 
        "the Power of the Count Five Test and F Test 
      for Different Sample Sizes.")
#Compare by standard values.
re_t.1 = (re_t[2,]-re_t[1,])/re_t[1,]
re_t.2 = (re_t[2,]-re_t[1,])/re_t[2,]
gap_t = round(data.frame(re_t.1,re_t.2),3)
colnames(gap_t) = c("gap/Count5test","gap/Ftest")
kable(gap_t, caption = 
        "the standard gap of Power of the Count Five Test 
      and F Test for Different Sample Sizes.")

## -----------------------------------------------------------------------------
options(warn = -1)
library(knitr)
library(boot)
library(dplyr)
data(aircondit)
data1 = aircondit$hours
MLE_Q1 = 1/mean(data1)
set.seed(1014)
MLE <- function(x,i) {
  1/mean(x[i])
}
obj = boot(data = data1, statistic = MLE, R = 1e4)
result_Q1 = round(c(original.MLE=obj$t0,bias=mean(obj$t)-obj$t0,
se=sd(obj$t)),3)
kable(data.frame(result_Q1),caption = "MLE of the Hazard Rate, Bias
and Standard Error of the Estimate")

## -----------------------------------------------------------------------------
set.seed(1014)
boot.mean <- function(x,i) mean(x[i])
de <- boot::boot(data = data1,statistic = boot.mean, R = 1e4)
ci <- boot.ci(de,type=c("norm","basic","perc","bca"))
ci

## -----------------------------------------------------------------------------
#compare the length of intervals.
len.basic = ci$basic[5]-ci$basic[4]
len.bca = ci$bca[5]-ci$bca[4]
len.perc = ci$percent[5]-ci$percent[4]
len.normal = ci$normal[3]-ci$normal[2]
result_Q2 = cbind(len.normal,len.basic,len.perc,len.bca)
kable(result_Q2,caption = "Interval Length of Different Methods")

## -----------------------------------------------------------------------------
hist(data1,breaks = 100,main = "Histogram of source data",xlab = "times in hours between failures")

## -----------------------------------------------------------------------------
set.seed(1014)
n_Q3 = 1e3
#Generate 1000 sample from a normal population with parameter mu=0,sigma=1.
sample = rnorm(n_Q3) 
de <- boot::boot(data = sample,statistic = boot.mean, R = 1e4)
ci <- boot.ci(de,type=c("norm","basic","perc"))
#show the confidence interval of method normal, basic and percentile.
ci

## -----------------------------------------------------------------------------
mu = 0 #true value of mean.
n = 5*1e1 #sample num.
m = 1e2 #num to calculate cover rate.
set.seed(1014)
ci.norm<-ci.basic<-ci.perc<-ci.bca<-matrix(NA,m,2)
for(i in 1:m){
sample = rnorm(n)
de = boot::boot(data=sample,statistic=boot.mean, R = 999)
ci = boot::boot.ci(de,type=c("norm","basic","perc"))
ci.norm[i,] = ci$norm[2:3]
ci.basic[i,] = ci$basic[4:5]
ci.perc[i,] = ci$percent[4:5]
}
cover = c(mean(ci.norm[,1]<=mu & ci.norm[,2]>=mu),mean(ci.basic[,1]<=mu & ci.basic[,2]>=mu),mean(ci.perc[,1]<=mu & ci.perc[,2]>=mu))
left = c(mean(ci.norm[,2]<=mu),mean(ci.basic[,2]<=mu),mean(ci.perc[,2]<=mu))
right = c(mean(ci.norm[,1]>=mu),mean(ci.basic[,1]>=mu),mean(ci.perc[,1]>=mu))
#codes below are to show the results in the form of table.
result_Q3 = matrix(c(cover,left,right),nrow = 3)
rownames(result_Q3) = c("Nomal","Basic","Percentile")
colnames(result_Q3) = c("in","left","right")
kable(result_Q3,caption = "Sample Mean Cover Rate in Different Position of 95% Confidence Intervals with Different Methods")

## -----------------------------------------------------------------------------
options(warn = -1)
library(knitr)
library(bootstrap)
data(scor)
data_Q1 = data.frame(scor)
colnames(data_Q1) = c("Mechanics","Vectors","Algebra","Analysis","Statistic")
n_Q1 = dim(data_Q1)[1]
#define function to compute MLE of Sigma and it's engenvalues.
MLE_theta <- function(data,n){
  #compute the MLE of Sigma.
  var_mle = (n-1)/n*var(data)
  #compute estimation of theta.
  lam = eigen(var_mle)$values
  theta.hat = lam[1]/sum(lam)
  return(theta.hat)
}

## -----------------------------------------------------------------------------
theta.hat = MLE_theta(data_Q1,n_Q1)
theta.jack = numeric(n_Q1)
for(i in 1:n_Q1){
theta.jack[i] <- MLE_theta(data_Q1[(1:n_Q1)[-i],],n_Q1-1)
}
bias.jack = (n_Q1 - 1) * (mean(theta.jack) - theta.hat)
se.jack = sqrt((n_Q1 - 1) * mean((theta.jack - theta.hat)^2))
#codes below are to show the result in the form of table.
result_Q1 = round(c(original_theta_hat = theta.hat,
                    bias.jack = bias.jack, se.jack = se.jack),3)
kable(data.frame(result_Q1), 
caption = "Estimation of Theta, Bias and Standard Error of the Estimate with Jackknife")

## -----------------------------------------------------------------------------
library(DAAG); attach(ironslag)
a <- seq(10, 40, .1) #sequence for plotting fits
#creat models.
#L1 <- lm(magnetic ~ chemical)
#L2 <- lm(magnetic ~ chemical + I(chemical^2))
#L3 <- lm(log(magnetic) ~ chemical)
#L4 <- lm(log(magnetic) ~ log(chemical))

n <- length(magnetic) #in DAAG ironslag
e1 <- e2 <- e3 <- e4 <- matrix(0,n,n)
# for n-fold cross validation
# fit models on leave-one-out samples
for (k in 1:n) {
  for (j in (1:n)[-k]) {
    y <- magnetic[-c(k,j)]
    x <- chemical[-c(k,j)]
    J1 <- lm(y ~ x)
    yhat1.1 <- J1$coef[1] + J1$coef[2] * chemical[k]
    yhat1.2 <- J1$coef[1] + J1$coef[2] * chemical[j]
    e1[k,j] <- ((magnetic[k] - yhat1.1)^2+(magnetic[j] - yhat1.2)^2)/2
    J2 <- lm(y ~ x + I(x^2))
    yhat2.1 <- J2$coef[1] + J2$coef[2] * chemical[k] +
    J2$coef[3] * chemical[k]^2
    yhat2.2 <- J2$coef[1] + J2$coef[2] * chemical[j] +
    J2$coef[3] * chemical[j]^2
    e2[k,j] <- ((magnetic[k] - yhat2.1)^2+(magnetic[j] - yhat2.2)^2)/2
    J3 <- lm(log(y) ~ x)
    logyhat3.1 <- J3$coef[1] + J3$coef[2] * chemical[k]
    yhat3.1 <- exp(logyhat3.1)
    logyhat3.2 <- J3$coef[1] + J3$coef[2] * chemical[j]
    yhat3.2 <- exp(logyhat3.2)
    e3[k,j] <- ((magnetic[k] - yhat3.1)^2+(magnetic[j] - yhat3.2)^2)/2
    J4 <- lm(log(y) ~ log(x))
    logyhat4.1 <- J4$coef[1] + J4$coef[2] * log(chemical[k])
    logyhat4.2 <- J4$coef[1] + J4$coef[2] * log(chemical[j])
    yhat4.1 <- exp(logyhat4.1)
    yhat4.2 <- exp(logyhat4.2)
    e4[k,j] <- ((magnetic[k] - yhat4.1)^2+(magnetic[j] - yhat4.2)^2)/2
  }
}
detach(ironslag)
#codes below are to show the result in the form of table.
result_Q2=c(sum(e1)/(n*(n-1)), 
            sum(e2)/(n*(n-1)), sum(e3)/(n*(n-1)), sum(e4)/(n*(n-1)))
table_Q2 = round(matrix(result_Q2),3)
rownames(table_Q2) = c("Model1","Model2","Model3","Model4")
colnames(table_Q2) = c("estimates for prediction error")
kable(table_Q2,
      caption = "Estimates for Prediction Error of Different Models")

## -----------------------------------------------------------------------------
result_Q2.one = c(19.55644, 17.85248, 18.44188, 20.45424)
differ = result_Q2 - result_Q2.one
comparison = round(matrix(c(result_Q2,result_Q2.one,differ),nrow = 4),3)
#codes below are to show the result in the form of table.
colnames(comparison) = c("leave-two-out","leave-one-out","difference")
rownames(comparison) = c("Model1","Model2","Model3","Model4")
kable(comparison,
      caption = "Estimates for Prediction Error of Different Models")

## -----------------------------------------------------------------------------
library(MASS)
set.seed(1021)
#generate data by multi-normal distribution.
mu = matrix(c(0,0))
sigma = matrix(c(4,1,1,1),nrow = 2)
data_Q3 = mvrnorm(n = 10,mu = mu,Sigma = sigma, tol = 1e-6)
x_Q3 = data_Q3[,1];y_Q3 = data_Q3[,2]
n_Q3 = length(x_Q3)
R = 999 #number of replicates
reps = numeric(R) #storage for replicates
t0 = cor(x_Q3, y_Q3, method = "spearman")
for (i in 1:R) {
#generate indices k for y.
k = sample(1:n_Q3, replace = FALSE)
x1 = x_Q3
y1 = y_Q3[k] #complement of x1
reps[i] = cor(x1, y1, method = "spearman")
}
p_per <- mean(abs(reps) >= abs(t0))

## -----------------------------------------------------------------------------
p_int = cor.test(x_Q3, y_Q3, method = "spearman", exact=FALSE)$p.value
result_Q3 = round(matrix(c(p_int,p_per)),3)
rownames(result_Q3) = c("pure cor.test","permutation test")
colnames(result_Q3) = c("p-value")
kable(result_Q3, table.envir="table*", 
      caption = "p-values of permutation test and pure Spearman rank correlation test")

## -----------------------------------------------------------------------------
options(warn = 1)
library(knitr)
library(VGAM)
rw.Metropolis <- function(x0, sigma, N){
  #Laplace distribution with parameter:loc=0,scale=1;
  #N: length of the chain;
  #x0: initial value
  #sigma:param of normal distribution.
  x = numeric(N)
  x[1] = x0
  u = runif(N)
  k = 0
  for (i in 2:N) {
    y = rnorm(1, x[i-1], sigma)
    if (u[i] <= (dlaplace(y,0,1) / dlaplace(x[i-1],0,1)))
      x[i] = y else {
      x[i] = x[i-1]
      k = k + 1
    }
  }
  return(list(x=x, k=k))
}

N = 2000
sigma = c(.05, .5, 2)
x0 = 25

#Four chains are generated for different variances sigma^2 of the proposal distribution.
set.seed(123)
rw1 = rw.Metropolis(x0, sigma[1], N)
rw2 = rw.Metropolis(x0, sigma[2], N)
rw3 = rw.Metropolis(x0, sigma[3], N)

#plot the four chains. 
#par(mfrow=c(2,2),oma = c(0, 0, 3, 0)) #display 4 graphs together
refline <- qlaplace(c(.025, .975))
rw <- cbind(rw1$x, rw2$x, rw3$x)
for (j in 1:3) {
plot(rw[,j], type="l",xlab=bquote(sigma == .(round(sigma[j],3))),ylab="X", ylim=range(rw[,j]))
abline(h=refline)
}
mtext("Random walk Metropolis chains with different variances", side = 3, line = 0, outer = T)
#par(mfrow=c(1,1)) #reset to default

## -----------------------------------------------------------------------------
#compare with the theoretical quantiles of the target distribution.
a = c(.05, seq(.1, .9, .1), .95)
Q = qlaplace(a)
rw = cbind(rw1$x, rw2$x, rw3$x)
#Discard the burn-in values in the first 500 rows of each chain.
mc = rw[501:N, ]
Qrw = apply(mc, 2, function(x) quantile(x, a))
result_Q1.2 = round(cbind(Q, Qrw), 3) 
colnames(result_Q1.2) = c("Q","rw1","rw2","rw3")
kable(result_Q1.2,caption = "Quantiles of Target
Distribution and Chains")

## -----------------------------------------------------------------------------
#number of candidate points rejected.
result_Q1.1 = data.frame((N-c(rw1$k, rw2$k, rw3$k))/N)
rownames(result_Q1.1) = c("rw1","rw2","rw3")
colnames(result_Q1.1) = "Acceptance Rate"
kable(result_Q1.1,caption = "Acceptance Rates of
Each Chain")

## -----------------------------------------------------------------------------
Gelman.Rubin <- function(psi){
  # psi[i,j] is the statistic psi(X[i,1:j])
  # for chain in i-th row of X
  psi <- as.matrix(psi)
  n <- ncol(psi)
  k <- nrow(psi)
  psi.means <- rowMeans(psi) #row means
  B <- n * var(psi.means) #between variance est.
  psi.w <- apply(psi, 1, "var") #within variances
  W <- mean(psi.w) #within est.
  v.hat <- W*(n-1)/n + (B/n) #upper variance est.
  r.hat <- v.hat / W #G-R statistic
  return(r.hat)
}

x0_GR = c(-10, -5, 5, 10)

k = 3 #number of chains to generate
n = 15000 #length of chains
b = 1000 #burn-in length
set.seed(123)
GR_R = numeric(3)
for (j in 1:3) {
  X = matrix(0, nrow=k, ncol=n)
  for (i in 1:k)
   X[i, ] = rw.Metropolis(x0_GR[i],sigma[j],n)$x
  
  #compute diagnostic statistics
  psi = t(apply(X, 1, cumsum))
  for (i in 1:nrow(psi))
  psi[i,] = psi[i,] / (1:ncol(psi))
  GR_R[j] = Gelman.Rubin(psi)
  
  #plot the sequence of R-hat statistics
  rhat = rep(0, n)
  for (j in (b+1):n)
  rhat[j] = Gelman.Rubin(psi[,1:j])
  plot(rhat[(b+1):n], type="l", xlab="", ylab="R",main="")
  abline(h=1.2, lty=2)
}
result_Q2.0 = round(data.frame(GR_R),3)
colnames(result_Q2.0) = c("Gelman Rubin R_hat")
rownames(result_Q2.0) =  c("rw1","rw2","rw3")
kable(result_Q2.0,caption = "Gelman Rubin R_hat for chains with diffenrent variance")

## -----------------------------------------------------------------------------
#initialize constants and parameters
N = 5000 #length of chain
burn = 1000 #burn-in length
X = matrix(0, N, 2) #the chain, a bivariate sample
rho = 0.9 #correlation
mu1 = 0
mu2 = 0
sigma1 = 1
sigma2 = 1
s1 = sqrt(1-rho^2)*sigma1
s2 = sqrt(1-rho^2)*sigma2
# generate the chain
X[1, ] = c(mu1, mu2) #initialize
Gibbs.ch <- function(N,sigma1,sigma2,mu1,mu2,s1,s2,X){
  for (i in 2:N) {
  x2 = X[i-1, 2]
  m1 = mu1 + rho * (x2 - mu2) * sigma1/sigma2
  X[i, 1] = rnorm(1, m1, s1)
  x1 = X[i, 1]
  m2 = mu2 + rho * (x1 - mu1) * sigma2/sigma1
  X[i, 2] = rnorm(1, m2, s2)
  }
  return(X)
}
set.seed(123)
X = Gibbs.ch(N,sigma1,sigma2,mu1,mu2,s1,s2,X)
b = burn + 1
x = X[b:N, ]
plot(x, main="", cex=.5, xlab=bquote(X[t]),ylab=bquote(Y[t]), ylim=range(x[,2]))

## -----------------------------------------------------------------------------
colnames(x) = c("X","Y")
model.lm = lm(Y~X,data = data.frame(x))
summary(model.lm)
plot(model.lm)

## -----------------------------------------------------------------------------

k=4
X = matrix(0, nrow=2*k, ncol=N)
#initial setting
X[c(1,1+k),1]=c(mu1+1, mu2+1)
X[c(2,2+k),1]=c(mu1-1, mu2-1)
X[c(3,3+k),1]=c(mu1-5, mu2-5)
X[c(4,4+k),1]=c(mu1+5, mu2+5)
set.seed(123)
for (i in 1:k)
  X[c(i,i+k), ] = t(Gibbs.ch(N,sigma1,sigma2,mu1,mu2,s1,s2,t(X[c(i,i+k),])))
  
#compute diagnostic statistics
X1 = X[1:4,]
X2 = X[5:8,]
psi.1 = t(apply(X1, 1, cumsum))
psi.2 = t(apply(X2, 1, cumsum))

for (i in 1:nrow(psi.1)){
  psi.1[i,] = psi.1[i,] / (1:ncol(psi.1))
  psi.2[i,] = psi.2[i,] / (1:ncol(psi.2))
}
result_Q2 = round(data.frame(c(Gelman.Rubin(psi.1),Gelman.Rubin(psi.2))),3)
rownames(result_Q2) = c("X_t","Y_t")
colnames(result_Q2) = c("Gelman.Rubin R_hat")
kable(result_Q2,caption = "Gelman-Rubin R_hat for X_t and Y_t")

#plot the sequence of R-hat statistics
rhat.1 = rhat.2 = rep(0, N)
for (j in b:N){
  rhat.1[j] = Gelman.Rubin(psi.1[,1:j])
  rhat.2[j] = Gelman.Rubin(psi.2[,1:j])
}
plot(rhat.1[b:N], type="l", xlab="", ylab="R",main=quote(X[t]))
abline(h=1.2, lty=2)
plot(rhat.2[b:N], type="l", xlab="", ylab="R",main=quote(Y[t]))
abline(h=1.2, lty=2)

## -----------------------------------------------------------------------------
options(warn = -1)
library(knitr)

## -----------------------------------------------------------------------------
n = 50
SE <- function(a,b,sa,sb){
  sqrt(a^2*sb+b^2*sa)
}

statis <- function(x,m,y){
  model1 = summary(lm(y~m+x))$coefficients
  model2 = summary(lm(m~x))$coefficients
  alpha = model2[2]
  se_alpha = model2[4]
  beta = model1[2]
  se_beta = model1[5]
  T = alpha*beta/SE(alpha,beta,se_alpha^2,se_beta^2)
  return(T)
}
#initial
p_per1 = p_per2 = p_per3 = 0

#testing
for (j in 1:500){
  X = rnorm(n,3,10)
  #a_M = 0, a_Y = 0, alpha = 0,beta=1,gamma = 1
  e_M = rnorm(n,0,1)
  e_Y = rnorm(n,0,1)
  M = e_M
  Y = X+M+e_Y

  #permutation of X and M
  R = 100
  reps = numeric(R)
  to = statis(X,M,Y)
  for (i in 1:R) {
    #generate indices k.
    k = sample(1:n, replace = FALSE)
    x1 = X[k]
    m1 = M
    y1 = Y
    model1 = summary(lm(Y~M+X))$coefficients
    model2 = summary(lm(m1~x1))$coefficients
    alpha = model2[2]
    se_alpha = model2[4]
    beta = model1[2]
    se_beta = model1[5]
    reps[i] = alpha*beta/SE(alpha,beta,se_alpha^2,se_beta^2)
  }
  p1 = mean(abs(reps) >= abs(to))
  if (p1<0.05) p_per1 = p_per1+1
}


## -----------------------------------------------------------------------------
for (j in 1:500){
  #a_M = 0, a_Y = 0, alpha = 1, beta=0,gamma = 1
  X = rnorm(n,3,10)
  #a_M = 0, a_Y = 0, alpha = 0,beta=1,gamma = 1
  e_M = rnorm(n,0,1)
  e_Y = rnorm(n,0,1)
  M = X+e_M
  Y = X+e_Y
  #permutation of Y and M
  reps = numeric(R)
  to = statis(X,M,Y)
  for (i in 1:R) {
    #generate indices k.
    k = sample(1:n, replace = FALSE)
    x1 = X
    y1 = Y[k]
    m1 = M
    model1 = summary(lm(y1~m1+x1))$coefficients
    model2 = summary(lm(M~X))$coefficients
    alpha = model2[2]
    se_alpha = model2[4]
    beta = model1[2]
    se_beta = model1[5]
    reps[i] = alpha*beta/SE(alpha,beta,se_alpha^2,se_beta^2)
  }
  p2 = mean(abs(reps) >= abs(to))
  if (p2<0.05) p_per2 = p_per2+1
} 


## -----------------------------------------------------------------------------
for (j in 1:500){
  #a_M = 0, a_Y = 0, alpha = 0, beta=0,gamma = 1
  X = rnorm(n,3,10)
  #a_M = 0, a_Y = 0, alpha = 0,beta=1,gamma = 1
  e_M = rnorm(n,0,1)
  e_Y = rnorm(n,0,1)
  M = e_M
  Y = X+e_Y
  #permutation of X, M and Y.
  reps = numeric(R)
  to = statis(X,M,Y)
  for (i in 1:R) {
    #generate indices k.
    k = sample(1:n, replace = FALSE)
    x1 = X
    m1 = M[k]
    y1 = Y
    model1 = summary(lm(y1~m1+x1))$coefficients
    model2 = summary(lm(m1~x1))$coefficients
    alpha = model2[2]
    se_alpha = model2[4]
    beta = model1[2]
    se_beta = model1[5]
    reps[i] = alpha*beta/SE(alpha,beta,se_alpha^2,se_beta^2)
  }
  p3 = mean(abs(reps) >= abs(to))
  if (p3<0.05) p_per3 = p_per3+1
}


## -----------------------------------------------------------------------------
result_Q1 = c(p_per1/500, p_per2/500, p_per3/500)
table_Q1 = round(matrix(result_Q1),5)
rownames(table_Q1) = c("alpha=0","beta=0","alpha=beta=0")
colnames(table_Q1) = "Type 1 error rate"
kable(table_Q1,caption = "Type 1 error rate for Different Occasions")

## -----------------------------------------------------------------------------
#define the function.
Q2func <- function(N,b1,b2,b3,f0){
  x1 = rpois(N,1) 
  x2 = rexp(N,1)
  x3 = sample(0:1,N,replace=TRUE)
  g = function(alpha){
  tmp = exp(-alpha-b1*x1-b2*x2-b3*x3)
  p = 1/(1+tmp)
  mean(p) - f0
  }
  solution = uniroot(g,c(-20,0))
  alpha = solution$root
  return(alpha)
}

## -----------------------------------------------------------------------------
f0 = c(0.1,0.01,0.001,0.0001)
Q2func1 <- function(f0){
  N = 1e6; b1 = 0; b2 = 1; b3 = -1
  Q2func(N,b1,b2,b3,f0)
}
set.seed(123)
result_Q2 = unlist(lapply(f0, Q2func1))
#codes below are to show the result in the form of table.
table_Q2 = round(matrix(c(f0,result_Q2),nrow = length(f0)),3)
colnames(table_Q2) = c("f0","alpha")
rownames(table_Q2) = NULL
kable(table_Q2,caption = "Table of Alpha vs F0")

## -----------------------------------------------------------------------------
plot(result_Q2, f0, pch = 16,xlab = "alpha",main = "f0 vs alpha")

## -----------------------------------------------------------------------------
options(warn = -1)
library(knitr)

## -----------------------------------------------------------------------------
utf8ToInt("one");utf8ToInt("2")

## -----------------------------------------------------------------------------
dim(c(1,2))

## -----------------------------------------------------------------------------
mat1 = matrix(c(1,2))
is.matrix(mat1);is.array(mat1)

## -----------------------------------------------------------------------------
data1 = data.frame(c(1,2))
attributes(data1)

## -----------------------------------------------------------------------------
data2 = data.frame(num = c(1,2),char = c("a","b"))
data2;as.matrix(data2)

## -----------------------------------------------------------------------------
#data frame with 0 rows
data3 = data.frame(num=integer())
data3;class(data3)
#data frame with 0 columns
data4 = data.frame(x=1)[,0,drop=FALSE]
dim(data4);class(data4)

## -----------------------------------------------------------------------------
# direct MLE
mle_Q4 <- function(U,V,lam){
  D = V-U
  C = sum(U)
  L = sum(D*exp(-D*lam)/(1-exp(-D*lam)))-C
  return(L)
}

E_step <- function(U,V,lam){
  s1 = -V*exp(-lam*V)-1/lam*exp(-lam*V)+U*exp(-lam*U)+1/lam*exp(-lam*U)
  s2 = 1-exp(-lam*V)-1+exp(-lam*U)
  s1/s2
}
M_step <- function(U,V,lam){
  x = E_step(U,V,lam)
  lam_new = length(U)/sum(x)
  if (abs(lam_new-lam)>0.0001) M_step(U,V,lam_new)
  else return(lam_new)
}
U = c(11,8,27,13,16,0,23,10,24,2)
V = c(12,9,28,14,17,1,24,11,25,3)
lam0 = 1 #initial value of EM.
mlE <- function(lam){
  mle_Q4(U,V,lam)
}
result_mle = uniroot(mlE,c(0,1))$root
result_EM = M_step(U,V,lam0)
result = round(matrix(c(result_mle,result_EM)),5)
rownames(result) = c("direct MLE","EM")
colnames(result)  = "lambda estimate"
kable(result,caption = "MLE of lambda using different algorithms")

## -----------------------------------------------------------------------------
options(warn = -1)

## -----------------------------------------------------------------------------
#function definition.
scale01 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  (x - rng[1]) / (rng[2] - rng[1])
}

data1 = data.frame(num1=c(1,5,10),num2=c(6,8,9))
apply(data1, 2, scale01)

## -----------------------------------------------------------------------------
#consider of the situation that there's any other type of data in some columns.
data2 = data.frame(num1=c(1,5,10),num2=c(6,8,9),name=c("A","B","C"))
#select columns that only contain numeric data.
num_cols = unlist(lapply(data2, is.numeric)) 
apply(data2[,num_cols], 2, scale01)

## -----------------------------------------------------------------------------
vapply(data1, sd, numeric(1))

## -----------------------------------------------------------------------------
#select columns that only contain numeric data.
data2
num_cols = vapply(data2, is.numeric, logical(1)) 
vapply(data2[,num_cols], sd, numeric(1))

## -----------------------------------------------------------------------------
library(Rcpp) 
library(microbenchmark)

## -----------------------------------------------------------------------------
# Define function "Gibbs_binorm"
# Can create source file in Rstudio
# sourceCpp("Gibbsbinorm.cpp")

## -----------------------------------------------------------------------------
#initialize constants and parameters
N = 5000 #length of chain
burn = 1000 #burn-in length
X = matrix(0, N, 2) #the chain, a bivariate sample
rho = 0.9 #correlation
mu1 = 0
mu2 = 0
sigma1 = 1
sigma2 = 1
s1 = sqrt(1-rho^2)*sigma1
s2 = sqrt(1-rho^2)*sigma2
# generate the chain
X[1, ] = c(mu1, mu2) #initialize
set.seed(1)
X = Gibbs_binorm(N,sigma1,sigma2,mu1,mu2,rho,s1,s2,X)
b = burn + 1
x = X[b:N, ]
plot(x, main="", cex=.5, xlab=bquote(X[t]),ylab=bquote(Y[t]), ylim=range(x[,2]))

## -----------------------------------------------------------------------------

Gibbs.ch <- function(N,sigma1,sigma2,mu1,mu2,rho,s1,s2,X){
  for (i in 2:N) {
  x2 = X[i-1, 2]
  m1 = mu1 + rho * (x2 - mu2) * sigma1/sigma2
  X[i, 1] = rnorm(1, m1, s1)
  x1 = X[i, 1]
  m2 = mu2 + rho * (x1 - mu1) * sigma2/sigma1
  X[i, 2] = rnorm(1, m2, s2)
  }
  return(X)
}
X1 = matrix(0, N, 2)
X1[1,] = c(mu1, mu2) #initialize
set.seed(2)
X1 = Gibbs.ch(N,sigma1,sigma2,mu1,mu2,rho,s1,s2,X1)
x1 = X1[b:N, ]
qqplot(x1[,1],x[,1],xlab = "pure R language-X",ylab = "Cpp-X")
qqplot(x1[,2],x[,2],xlab = "pure R language-Y",ylab = "Cpp-Y")

## -----------------------------------------------------------------------------
ts <- microbenchmark(GibbsR=Gibbs.ch(N,sigma1,sigma2,mu1,mu2,rho,s1,s2,X),Gibbscpp=Gibbs_binorm(N,sigma1,sigma2,mu1,mu2,rho,s1,s2,X))
summary(ts)[,c(1,3,5,6)]


