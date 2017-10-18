library(testthat)
library(ESSpackage)
library(MVN)
library(MASS)
library(mvtnorm)
library(ICS)
library(ggplot2)

####################################################################
############# COMPARISON FOR 2-DIM GAUSSIAN REGRESSION #############
####################################################################

test_that("ESSpackage",{

  ## log likelihood function
  gr = function(yn,std){
    function(f){
      sum(dnorm(yn,f,std,log = TRUE))
    }
  }

  # Prior setup
  x_toy = matrix(runif(2),1,2)
  sigma_toy = cov_mat(1,1,x_toy,2)
  f_toy_new = mvrnorm(n = 1, mu = rep(0,dim(sigma_toy)[1]), sigma_toy)

  #Likelihood setup
  std_gr = 0.3
  sd_lik <- matrix(c(0.09,0,0,0.09),2,2)
  observation_toy_new = c(rnorm(1,f_toy_new[1],sd = 0.3),rnorm(1,f_toy_new[2],sd = 0.3))

  # Posterior setup
  mu_post <- solve(solve(sigma_toy)+solve(sd_lik))%*%solve(sd_lik)%*%(observation_toy_new)
  sd_post <- solve(solve(sigma_toy)+solve(sd_lik))

  # 4 models: Theoretical, ESS, NealMH, adaptiveMH
  fpost_toy <- mvrnorm(100000, mu_post, sd_post)
  r_toy_new = ess(f=f_toy_new,sigma=sigma_toy,llk=gr(yn = observation_toy_new,std = std_gr),n=1000)
  neal_toy <- NMH(f=f_toy_new,sigma=sigma_toy,llk=gr(yn = observation_toy_new,std = std_gr),1000,0.2)
  AdaptMH_toy <- AdaptMH(f=f_toy_new,sigma=sigma_toy,llk=gr(yn = observation_toy_new,std = std_gr),n=5000,N=2,cov_post=sd_post)



  # Contour Plot of the 3 models
  df <- data.frame(fpost_toy)
  df_1 <- data.frame(r_toy_new)
  df_2 <- data.frame(neal_toy$fs)
  df_3 <- data.frame(AdaptMH_toy)
  library(ggplot2)
  p <- ggplot() + geom_density2d(aes(x = X1, y = X2, colour="green"), data = df) +
    geom_density2d(aes(x = X1, y = X2, colour="blue"), data = df_1) +
    geom_density2d(aes(x = X1, y = X2, col="red"), data = df_2) +
    geom_density2d(aes(x = X1, y = X2, col="yellow"), data = df_3) +
    scale_colour_manual(values=c("green", "blue", "red", "yellow"),labels = c("Theoretical", "ESS", "NealMH", "AdaptiveMH"), name="Sampling methods")
  print(p + ggtitle("Comparison of contour plots"))

})
