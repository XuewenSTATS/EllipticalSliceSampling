# Elliptical Slice Sampling
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'


## elliptical slice sampler
ess = function(f,sigma,llk,n){
  fs = matrix(0,n,length(f))
  fs[1,] = f
  for(j in 1:(n-1)) {
    nu = mvrnorm(n = 1, mu = rep(0,dim(sigma)[1]), sigma)
    logy = log(runif(1)) + llk(f)
    theta = runif(1,0,2 * pi)
    min = theta - 2 * pi
    max = theta
    f1 = f * cos(theta) + nu * sin(theta)
    ## loop
    while(llk(f1) < logy){
      if(theta < 0){
        min = theta
      }
      else{
        max = theta
      }
      theta = runif(1,min,max)
      f1 = f * cos(theta) + nu * sin(theta)
    }
    fs[j+1,] = f1
    f = f1
    print(j)
  }
  return(fs = fs)
}


## Covariance matrix for f
## sig_var is the unit signal variance, default = 1
## l is the lengthscale, default = 1
## N is the number of columns, which is 200 in default
cov_mat = function(sig_var ,l,x,N){
  mat = matrix(0,N,N)
  for (i in 1:dim(x)[2]){
    for(j in 1:i){
      mat[i,j] = sig_var * exp(-0.5*(t(x[,i]-x[,j])%*%(x[,i]-x[,j]))/l^2)
      mat[j,i] = mat[i,j]
    }
  }
  return(mat)
}


## one possible log likelihood
make.mvn <- function(mean, vcv) {
  logdet <- as.numeric(determinant(vcv, TRUE)$modulus)
  tmp <- length(mean) * log(2 * pi) + logdet
  vcv.i <- solve(vcv)

  function(x) {
    dx <- x - mean
    -(tmp + rowSums((dx %*% vcv.i) * dx))/2
  }
}

## library(mvnormtest)
## mshapiro.test(t(ess(f,sigma,llk = make.mvn(c(0,0),sigma),n=100)))
## not work for univariate gaussian (checked by shapiro test)

## check it on Gaussian regression

## try 1 dimension data
library(MASS)
## parameters: (gr:gaussian regression)
N_gr = 200  ## sample size
var_gr = 0.09  ## noise variance
std_gr = 0.3  ## standard deviation
l_gr = 1
sig_var_gr = 1
## inputs 
x1 = matrix(runif(N_gr),nrow = 1,ncol = N_gr)
sigma1 = cov_mat(sig_var = sig_var_gr,l = l_gr,x = x1,N=N_gr)
f1 = mvrnorm(n = 1, mu = rep(0,dim(sigma1)[1]), sigma1)
## generate observations
observation1 = numeric(N_gr)
for(j in 1:N_gr){
  observation1[j] = rnorm(1,f1[j],sd = std_gr)
}
## log likelihood function
gr = function(yn,std){
  function(f){
    sum(dnorm(yn,f,std,log = TRUE))
  }
}

## toy example
## before we implement the code on 1 dimension data, we first set N = 2, D = 1
x_toy = matrix(runif(2),1,2)
sigma_toy = cov_mat(1,1,x_toy,2)
f_toy = mvrnorm(n = 1, mu = rep(0,dim(sigma_toy)[1]), sigma_toy)
observation_toy = c(rnorm(1,f_toy[1],sd = 0.3),rnorm(1,f_toy[2],sd = 0.3))
r_toy = ess(f=f_toy,sigma=sigma_toy,llk=gr(yn = observation_toy,std = std_gr),n=100000)
library(MVN)
pdf("qqplot.pdf",6,4)
uniPlot(r_toy[-(1:10000),],type = "qqplot")
dev.off()
pdf("histogram.pdf",6,4)
uniPlot(r_toy[-(1:10000),],type = "histogram")
dev.off()
results1 = mardiaTest(r_toy[-(1:80000),],qqplot = FALSE)
pdf("3dplot.pdf",6,4)
par(mfrow = c(1,2))
mvnPlot(results1, type = "persp", default = TRUE)
mvnPlot(results1, type = "contour", default = TRUE,xlab="f_1",ylab="f_2")
dev.off()

## start with 1 dimensional example
r1 = ess(f=f1,sigma=sigma1,llk=gr(yn = observation1,std = std_gr),n=1000)
system.time(ess(f=f1,sigma=sigma,llk=gr(yn = observation,std = std_gr),n=1000))  ## 20.148
llk=gr(yn = observation1,std = std_gr)
loglikr1 = apply(r1,1,llk)
## thin it (every 3 iterations)
loglikr1_thined = numeric(333)
iterations = numeric(333)
for(k in 1:333){
  loglikr1_thined[k] = loglikr1[3*k]
  iterations[k] = 3*k
}

pdf("figure1.pdf",5,4)
plot(iterations,loglikr1_thined,type="l",col="blue",xlab = "# iterations",ylab = "",main = "Elliptical slice sampling")
dev.off()

## convergence diagnstic
r1_conv = ess(f=f1,sigma=sigma1,llk=gr(yn = observation1,std = std_gr),n=100000)

## 1. trace plot
plot(r1_conv[1:50000,1],type="l",col= "blue",ylab= "f_1",xlab = "iterations")
plot(apply(r1_conv[1:50000,],1,llk),type="l",col= "blue",ylab= "log likelihood",xlab = "iterations")
## (?? I cant observe burn-in from traceplot)
## 2. effective sample size
library(coda)
effectiveSize(apply(r1_conv[-(1:10000),],1,llk))  ## 4077
effectiveSize(r1_conv[-(1:10000),1]) ## 3177
## 3. check normality
install.packages("MVN")
library(MVN)
uniPlot(r1_conv[-(1:10000),1:4],type = "qqplot")
uniPlot(r1_conv[-(1:10000),1:4],type = "histogram")
results1 = hzTest(r1_conv[90000:100000,1:2],qqplot = FALSE)
mvnPlot(results1, type = "persp", default = TRUE)
mvnPlot(results1, type = "contour", default = TRUE)
library(mvnormtest)
mshapiro.test(r1_conv[90000:100000,])
## 100 runs of 100000 iterations with 10000 burn-in

t1 = numeric(100)
ES = numeric(100)
for(j in 1:100){
  t1[j] = system.time({R1 = ess(f=f1,sigma=sigma1,llk=gr(yn = observation1,std = std_gr),n=100000)})[3]
  loglikR1 = apply(R1,1,llk)
  ES[j] = effectiveSize(loglikR1[-(1:10000)])
}

## try 10 dimensions
x10 = matrix(runif(10 * N_gr),nrow = 10,ncol = N_gr)
sigma10 = cov_mat(sig_var = sig_var_gr,l = l_gr,x = x10,N=N_gr)
f10 = mvrnorm(n = 1, mu = rep(0,dim(sigma10)[1]), sigma10)
observation10 = numeric(N_gr)
for(j in 1:N_gr){
  observation10[j] = rnorm(1,f10[j],sd = std_gr)
}

r10 = ess(f=f10,sigma=sigma10,llk=gr(yn = observation10,std = std_gr),n=10000)
system.time(ess(f=f10,sigma=sigma10,llk=gr(yn = observation10,std = std_gr),n=10000))
## 244.375
llk10=gr(yn = observation10,std = stg_gr)
loglikr10 = apply(r10,1,llk10)
## thin it (every 78 iterations)
loglikr10_thined = numeric(128)
iterations = numeric(128)
for(k in 1:128){
  loglikr10_thined[k] = loglikr10[78*k]
  iterations[k] = 78*k
}

pdf("figure10.pdf",5,4)
plot(iterations,loglikr10_thined,type="l",col="blue",ylim = c(-150,-30),xlab = "# iterations",ylab = "",main = "Elliptical slice sampling")
dev.off()

## convergence diagnostic
r10_conv = ess(f=f10,sigma=sigma10,llk=gr(yn = observation10,std = std_gr),n=100000)

## 100 runs of 100000 iterations with 10000 burn-in
library(coda)
t10 = numeric(100)
ES10 = numeric(100)
for(j in 1:100){
  t10[j] = system.time({R10 = ess(f=f10,sigma=sigma10,llk=gr(yn = observation10,std = std_gr),n=100000)})[3]
  loglikR10 = apply(R1,1,llk10)
  ES10[j] = effectiveSize(loglikR1[-(1:10000)])
}



## log Gaussian cox process
mining = read.table("mining.dat")
str(mining)
library(pracma)
## parameters:
bin_width = 50
N_lgcp = 811 ## number of bins
sig_var_lgcp = 1
l_lgcp = 13516

get_mine_data = function(bin_width){
intervals = as.matrix(mining) ## 190 * 1 matrix
num_days = 40550
num_events = 191
event_days = c(1,cumsum(intervals)+1)
event_days[length(event_days)]
edges = c(seq(1,num_days,bin_width),num_days+1)
bin_counts = histc(event_days,edges)$cnt
if(bin_counts[length(bin_counts)] == 0){
  bin_counts = bin_counts[-length(bin_counts)]
}
xx = (edges[1:(length(edges)-1)] + (edges[2:length(edges)]-1))/2
yy = bin_counts
return(list(xx = xx,yy = yy))
}

llk_lgcp = function(yn,m){
  function(f){
  sum(dpois(yn,exp(f+m),log = TRUE))
  }
}
observation_lgcp = get_mine_data(bin_width = 50)$yy
## simulate f (length = 811)
x_lgcp = matrix(runif(N_lgcp),nrow = 1,ncol = N_lgcp)
sigma_lgcp = cov_mat(sig_var = sig_var_lgcp,l = l_lgcp,x_lgcp,N=N_lgcp)
f_lgcp = mvrnorm(n = 1, mu = rep(0,dim(sigma_lgcp)[1]), sigma_lgcp)
## parameters
m = log(num_events/N_lgcp)
llk = llk_lgcp(yn = observation_lgcp,m=m)
llk(f)
## convergence diagnostic
logGauss = ess(f=f_lgcp,sigma=sigma_lgcp,llk=llk_lgcp(yn = observation_lgcp,m=m),n=10000)
loglikrmine = apply(logGauss,1,llk)
pdf("figureminedata.pdf",5,4)
plot(iterations,loglikrmine_thined,type="l",col="blue",xlab = "# iterations",ylab = "",main = "Elliptical slice sampling")
dev.off()
plot(logGauss, type='l')

 




