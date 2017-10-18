## Packages
library(ESSpackage)
library(MASS)
library(MVN)
library(mvtnorm)
library(parallel)
library(foreach)
library(doSNOW)
library(coda)
library(pracma)

## In our simulation, we realised that the code was extremely slow for input with large dimensions.
## Hence we parallel the function in the package and output log likelihood only.

ess_llk = function(f,sigma,llk,n){

  foreach(j = 1:(n-1), .combine = "c",.packages = "MASS") %dopar% {
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
    f = f1
    llk(f)
  }
}

NMH_llk = function(f,sigma,llk,n,stepsize){
  foreach(j = 1:(n-1),.combine = "c",.packages = "MASS") %dopar% {
    nu = mvrnorm(n = 1, mu = rep(0,dim(sigma)[1]), sigma)
    logy = log(runif(1)) + llk(f)
    f1 = sqrt(1-stepsize^2) * f + stepsize * nu
    if(llk(f1) > logy) {
      f1 = f1
    }
    else{
      f1 = f
    }
    f = f1
    llk(f)
  }
}

#########################################################
################ Gaussian regression ####################
#########################################################

## log likelihood function
gr = function(yn,std){
  function(f){
    sum(dnorm(yn,f,std,log = TRUE))
  }
}

## toy example
## before we validate the code on 1 dimension data, we first set N = 2, D = 1
x_toy = matrix(runif(2),1,2)
sigma_toy = cov_mat(1,1,x_toy,2)
f_toy = mvrnorm(n = 1, mu = rep(0,dim(sigma_toy)[1]), sigma_toy)
observation_toy = c(rnorm(1,f_toy[1],sd = 0.3),rnorm(1,f_toy[2],sd = 0.3))
r_toy = ess(f=f_toy,sigma=sigma_toy,llk=gr(yn = observation_toy,std = 0.3),n=100000)
# convergence diagnostic
plot(r_toy[,1],type="l",col = "blue",ylab= "f_1",xlab = "iterations")
plot(r_toy[,2],type="l",col = "blue",ylab= "f_2",xlab = "iterations")
llk=gr(yn = observation_toy,std = 0.3)
plot(apply(r1_conv,1,llk),type="l",col = "blue",ylab= "f_2",xlab = "iterations")
acf(apply(r1_conv,1,llk))
effectiveSize(apply(r1_conv,1,llk))
# normality check
pdf("qqplot.pdf",6,4)
uniPlot(r_toy[-(1:10000),],type = "qqplot")
dev.off()
results1 = hzTest(r_toy[80000:100000,],qqplot = FALSE)
pdf("3dplot.pdf",6,4)
mvnPlot(results1, type = "persp", default = TRUE)
dev.off()

## 1 dimensional synthetic data
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

## before we run NMH_llk, we try on smaller number of iterations, 100000 on NMH to find a proper stepsize
## such that the Markov chain converge faster and normality is reached
## stepsize = 0.2 is found to be better
library(parallel)
cl = makeCluster(4)
library(foreach)
library(doSNOW)
registerDoSNOW(cl)
llk=gr(yn = observation1,std = std_gr)
clusterExport(cl,list("observation1","llk","std_gr"))
system.time({r1llk = ess_llk(f=f1,sigma=sigma1,llk=llk,n=1000000)})
system.time({r1llk1 = NMH_llk(f=f1,sigma=sigma1,llk=llk,n=1000000,stepsize = 0.2)})
stopCluster(cl)

effectiveSize(r1llk[-(1:100000)])
effectiveSize(r1llk1[-(1:100000)])


## 10 dimensional synthetic data
x10 = matrix(runif(10 * N_gr),nrow = 10,ncol = N_gr)
sigma10 = cov_mat(sig_var = sig_var_gr,l = l_gr,x = x10,N=N_gr)
f10 = mvrnorm(n = 1, mu = rep(0,dim(sigma10)[1]), sigma10)
observation10 = numeric(N_gr)
for(j in 1:N_gr){
  observation10[j] = rnorm(1,f10[j],sd = std_gr)
}

cl = makeCluster(4)
registerDoSNOW(cl)
clusterExport(cl,list("observation10","llk","std_gr"))
system.time({r10llk = ess_llk(f=f10,sigma=sigma10,llk=llk,n=1000000)})
system.time({r10llk1 = NMH_llk(f=f10,sigma=sigma10,llk=llk,n=1000000,stepsize = 0.2)})
stopCluster(cl)

effectiveSize(r10llk[-(1:100000)])
effectiveSize(r10llk1[-(1:100000)])

## log Gaussian cox process
mining = read.table("mining.dat")
## parameters:
bin_width = 50
N_lgcp = 811 ## number of bins
sig_var_lgcp = 1
l_lgcp = 13516
num_days = 40550
num_events = 191

get_mine_data = function(bin_width,num_days,num_events){
  intervals = as.matrix(mining) ## 190 * 1 matrix
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
observation_lgcp = get_mine_data(bin_width = bin_width,num_days = 40550,num_events = 191)$yy
## simulate f (length = 811)
x_lgcp = matrix(runif(N_lgcp),nrow = 1,ncol = N_lgcp)
sigma_lgcp = cov_mat(sig_var = sig_var_lgcp,l = l_lgcp,x_lgcp,N=N_lgcp)
f_lgcp = mvrnorm(n = 1, mu = rep(0,dim(sigma_lgcp)[1]), sigma_lgcp)
## parameters
m = log(num_events/N_lgcp)
llk = llk_lgcp(yn = observation_lgcp,m=m)
llk(f)

# system.time({logGauss = ess(f=f_lgcp,sigma=sigma_lgcp,llk=llk_lgcp(yn = observation_lgcp,m=m),n=10000)})
## we found this is very slow because of the high dimension of f
## so we modify the input by setting
bin_width = 400
N_lgcp = 102
observation_lgcp = get_mine_data(bin_width = bin_width,num_days = 40550,num_events = 191)$yy
## simulate f (length = 102)
x_lgcp = matrix(runif(N_lgcp),nrow = 1,ncol = N_lgcp)
sigma_lgcp = cov_mat(sig_var = sig_var_lgcp,l = l_lgcp,x_lgcp,N=N_lgcp)
f_lgcp = mvrnorm(n = 1, mu = rep(0,dim(sigma_lgcp)[1]), sigma_lgcp)
m = log(num_events/N_lgcp)

cl = makeCluster(4)
registerDoSNOW(cl)
clusterExport(cl,list("observation_lgcp","llk","sigma_lgcp"))
system.time({logGauss = ess_llk(f=f_lgcp,sigma=sigma_lgcp,llk=llk_lgcp(yn = observation_lgcp,m=m),n=1000000)})
system.time({logGauss2= NMH_llk(f_lgcp,sigma=sigma_lgcp,llk=llk_lgcp(yn = observation_lgcp,m=m),n=1000000,stepsize = 0.2)})
stopCluster(cl)
acf(logGauss)
# system.time({logGauss3 = AdaptMH_loglik(f_lgcp,sigma_lgcp,llk=llk_lgcp(yn = observation_lgcp,m=m),n=100,N=101,beta=0.05)})
# this failed because the empirical estimate of covariance matrix contains many NA entries

effectiveSize(logGauss[-(1:100000)])
effectiveSize(logGauss2[-(1:100000)])


