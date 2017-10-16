#' Elliptical Slice Sampling
#'
#' This is an example function named 'ess'
#' which implement Elliptical Slice Sampling algorithm.
#'
#' Elliptical Slice Sampling
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
#' Covariance Matrix Elliptical Slice Sampling
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













