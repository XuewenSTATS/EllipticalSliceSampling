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

## ess
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
  }
  return(fs = fs)
}

## library(mvnormtest)
## mshapiro.test(t(ess(f,sigma,llk = make.mvn(c(0,0),sigma),n=100)))

## not work for univariate gaussian (checked by shapiro test)

## check it on Gaussian regression
## Covariance matrix for f
            ## sig_var is the unit signal variance, default = 1
            ## l is the lengthscale, default = 1
            ## N is the number of columns, which is 200 in default
cov_mat = function(sig_var = 1,l = 1,x,N){
  mat = matrix(0,N,N)
  for (i in 1:dim(x)[2]){
    for(j in 1:i){
    mat[i,j] = sig_var * exp(-0.5*(t(x[,i]-x[,j])%*%(x[,i]-x[,j]))/l^2)
    mat[j,i] = mat[i,j]
    }
  }
  return(mat)
}

## simulate synthetic data x
## try 1 dimension data
library(MASS)
N = 200
set.seed(1)
x = matrix(runif(200),nrow = 1,ncol = N)
sigma = cov_mat(sig_var = 1,l = 1,x,N=200)
f = mvrnorm(n = 1, mu = rep(0,dim(sigma)[1]), sigma)
## Gaussian regression log-likelihood function
## parameters: N = 200
##             noise variance: var = 0.09, standard deviation = 0.3
observation = numeric(200)
for(j in 1:200){
  observation[j] = rnorm(1,f[j],sd = 0.3)
}
gr = function(yn,std){
  function(f){
    sum(dnorm(yn,f,std,log = TRUE))
  }
}
r1 = ess(f=f,sigma=sigma,llk=gr(yn = observation,std = 0.3),n=1000)
system.time(ess(f=f,sigma=sigma,llk=gr(yn = observation,std = 0.3),n=1000))
## 20.148
llk=gr(yn = observation,std = 0.3)
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

## try 10 dimensions
set.seed(10)
x10 = matrix(runif(2000),nrow = 10,ncol = N)
sigma10 = cov_mat(sig_var = 1,l = 1,x10,N=200)
f10 = mvrnorm(n = 1, mu = rep(0,dim(sigma10)[1]), sigma10)
observation10 = numeric(200)
for(j in 1:200){
  observation10[j] = rnorm(1,f10[j],sd = 0.3)
}
gr = function(yn,std){
  function(f){
    sum(dnorm(yn,f,std,log = TRUE))
  }
}
r10 = ess(f=f10,sigma=sigma10,llk=gr(yn = observation10,std = 0.3),n=10000)
system.time(ess(f=f10,sigma=sigma10,llk=gr(yn = observation10,std = 0.3),n=10000))
## 244.375
llk=gr(yn = observation10,std = 0.3)
loglikr10 = apply(r10,1,llk)
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

## log Gaussian cox process
mining = read.table("mining.dat")
str(mining)
mining = as.matrix(mining) ## 190 * 1 matrix
llk_lgcp = function(yn,m){
  sum(dpois(yn,exp(f+m),log = TRUE)) + log(m)
}
observation_lgcp = mining 
## parameters
m = log(191/811)
## Effective Sample Size
## read mining data

## 100000 iterations
# 1-dim
library(coda)
t1 = numeric(100)
ES = numeric(100)
llk=gr(yn = observation,std = 0.3)
for(j in 1:100){
t1[j] = system.time({R1 = ess(f=f,sigma=sigma,llk=gr(yn = observation,std = 0.3),n=100000)})[3]
loglikR1 = apply(R1,1,llk)
ES[j] = effectiveSize(loglikR1[-(1:10000)])
}




