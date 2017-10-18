#' Elliptical Slice Sampling
#'
#' This is an example function named 'ess'
#' which implement Elliptical Slice Sampling algorithm.
#'
#'
#'
#' @param f initial state
#' @param sigma covariance matrix of f
#' @param llk log likelihood function
#' @param n number of iterations
#' @return an n x N matrix where N is the dimension of f
#' @export f

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

#' Covariance Matrix of Elliptical Slice Sampling
#'
#' This is an example function named 'cov_mat'
#' which gives the covariance matrix of f, the vector of latent variables.
#'
#'
#' @param sig_var signal variance
#' @param l lengthscale parameter
#' @param x feature dataset
#' @param N size of data x
#' @return covariance matrix of f
#' @export mat
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


#' Neal's Metropolis Hastings (1999)
#'
#' This is a function named 'NMH'
#' which implement the Metropolis Hastings algorithm adapted by Neal in 1999.
#'
#'
#'
#' @param f initial state
#' @param sigma covariance matrix of f
#' @param llk log likelihood function
#' @param n number of iterations
#' @param stepsize this is the fixed number considered in proposal
#' @return This function returns a list which contains an n x N matrix where N is the dimension of f and the corresponding log likelihood.
#' @export f

NMH = function(f,sigma,llk,n,stepsize){
  fs = matrix(0,n,length(f))
  fs[1,] = f
  loglik = numeric(n)
  loglik[1] = llk(f)
  for(j in 1:(n-1)){
    nu = mvrnorm(n = 1, mu = rep(0,dim(sigma)[1]), sigma)
    logy = log(runif(1)) + llk(f)
    f1 = sqrt(1-stepsize^2) * f + stepsize * nu
    if(llk(f1) > logy) {
      fs[j+1,] = f1
      loglik[j+1] = llk(f1)
    }
    else{
      fs[j+1,] = fs[j,]
      loglik[j+1] = loglik[j]
    }
    f = fs[j+1,]
  }
  return(list(fs = fs, loglik = loglik))
}


#' Adaptive Metropolis Hastings
#'
#' This is a function named 'AdaptMH'
#' which implement the Metropolis Hastings algorithm proposed by Gareth and Jeffrey (2012). However, instead of
#' using proposals involving empirical estimate of covariance matrix, we input the known covariance of posterior distribution.
#'
#' We need to load package "mvtnorm" before using this function
#'
#' @param f initial state
#' @param sigma covariance matrix of f
#' @param N dimension of f
#' @param llk log likelihood function
#' @param n number of iterations
#' @param cov_post this is the theoretical covariance matrix of the posterior
#' @return an n x N matrix where N is the dimension of f
#' @export f

AdaptMH = function(f,sigma,llk,n,N,cov_post){
  fs = matrix(0,n,length(f))
  fs[1,] = f
  for(j in 1:(n-1)){
    f1 = mvrnorm(n = 1, mu = f, (2.38^2) * cov_post / N)
    mhr = dmvnorm(x = f1,mean = rep(0,N),sigma,log = TRUE) + llk(f1) - dmvnorm(x = f,mean = rep(0,N),sigma,log=TRUE) - llk(f)
    if(log(runif(1)) < mhr){
      fs[j+1,] = f1
    }
    else{
      fs[j+1,] = f
    }
    f = fs[j+1,]
    print(j)
  }
  return(fs = fs)
}


#' Adaptive Metropolis Hastings using proposals involving empirical estimate of covariance
#'
#' This is a function named 'AdaptMH_full'
#' which implement the Metropolis Hastings algorithm proposed by Gareth and Jeffrey (2012).
#'
#' We need to load package "mvtnorm" before using this function
#'
#' @param f initial state
#' @param sigma covariance matrix of f
#' @param N dimension of f
#' @param llk log likelihood function
#' @param n number of iterations
#' @param N dimension of f
#' @param beta fixed constant used in proposal, the default value is 0.05
#' @return This function returns a list which contains an n x N matrix where N is the dimension of f and the corresponding log likelihood.
#' @export f

AdaptMH_full = function(f,sigma,llk,n,N,beta=0.05){

  fs = matrix(0,n,length(f))
  fs[1,] = f

  loglik = numeric(n)
  loglik[1] = llk(f)
  if((n-1) <= 2 * N){
    for(j in 1:(n-1)){
      f1 = mvrnorm(n = 1, mu = f, (0.1^2) * diag(N) / N)
      mhr = log(dmvnorm(x = f1,mean = rep(0,N),sigma) * dmvnorm(x = f,mean = rep(0,N),sigma)) + llk(f1) - llk(f)
      if(log(runif(1)) < mhr){
        loglik[j+1] = llk(f1)
        fs[j+1,] = f1
      }
      else{
        loglik[j+1] = loglik[j]
        fs[j+1,] = f
      }
      f = fs[j+1,]
      print(j)
    }
  }
  else{
    for(j in 1:(2*N)){
      f1 = mvrnorm(n = 1, mu = f, (0.1^2) * diag(N) / N)
      mhr = log(dmvnorm(x = f1,mean = rep(0,N),sigma) * dmvnorm(x = f,mean = rep(0,N),sigma)) + llk(f1) - llk(f)
      if(log(runif(1)) < mhr){
        loglik[j+1] = llk(f1)
        fs[j+1,] = f1
      }
      else{
        loglik[j+1] = loglik[j]
        fs[j+1,] = f
      }
      f = fs[j+1,]
      print(j)
    }
    for(j in (2*N+1):(n-1)){

      v1 = mvrnorm(n = 1, mu = f, (2.38^2) * cor(fs[1:j-1,]) / N)
      v2 = mvrnorm(n = 1, mu = f, (0.1^2) * diag(N) / N)
      f1 = (1-beta) * v1 + beta * v2
      mhr = log(dmvnorm(x = f1,mean = rep(0,N),sigma) * dmvnorm(x = f,mean = rep(0,N),sigma)) + llk(f1) - llk(f)
      if(log(runif(1)) < mhr){
        loglik[j+1] = llk(f1)
        fs[j+1,] = f1
      }
      else{
        loglik[j+1] = loglik[j]
        fs[j+1,] = f
      }
      f = fs[j+1,]
      print(j)
    }
  }
  return(list(fs = fs, loglik = loglik))
}
