#' Panos Toulis (panos.toulis@chicagobooth.edu)
#' Invariant inference via residual randomization.
#' 
#' June, 2022
#' Reproduces Section 6.1 of paper, Table 6
#' 
rm(list=ls())
library(sandwich)
library(tidyverse)
library(expm)

RAND_REPS = 3000
get_R = function(x, z) {
  k = ifelse(is.null(ncol(z)), 1, ncol(z))
  R = matrix(c(0, 1, rep(0, k)), nrow=1)
  # R = c(0, 1, 0)
  return(R)
}
  
#' Test statistic
tstat = function(y, x, z) {
  R = get_R(x,z)
  fit = lm(y ~ x+z)
  V2b = vcovHC(fit)
  bhat = coef(fit)
  bR = R %*% bhat
  as.numeric(t(bR) %*% solve(R %*% V2b %*% t(R)) %*% (bR)) 
}

#' RRI using a linear test statistic.
#' @param use_restricted TRUE to use restricted OLS residuals.
RRI_lin = function(y, x, z, use_restricted=TRUE) {
  
  fit = lm(y ~ x+z)
  bhat = coef(fit)[2]
  
  X = cbind(1, x, z)
  Q = solve(t(X) %*% X) %*% t(X)
  R = get_R(x, z)
  aX = as.numeric(R %*% Q)
  # A_sqrt = sqrtm(A)
  
  tn = function(e) {
    #norm(A_sqrt %*% e, "2")^2
    sum(aX * e)
  }
  
  # tobs
  f0 = lm(y ~ z) # under H0
  e0 = f0$residuals
  tobs = bhat
  n = length(e0)
  
  if(!use_restricted) {
    e0 = fit$residuals
  }
  
  # randomization distribution
  tvals = replicate(RAND_REPS, {
    s = sign(rnorm(n))
    tn(s*sample(e0))  
    # tn(e*s)
  })
  
  p = mean(tvals >= tobs)
  2*min(p, 1-p) <= 0.05
}

#' Residual randomization using Quadratic test statistic
RRI_quad = function(y, x, z, use_restricted=TRUE) {
  
  R = get_R(x,z)
  fit = lm(y ~ x+z)
  bhat = coef(fit)
  V2b = vcovHC(fit)
  X = cbind(1, x, z)
  Q = solve(t(X) %*% X) %*% t(X)
  
  A = t(Q) %*% t(R) %*% solve(R%*% V2b %*% t(R)) %*% R %*% Q
  # A_sqrt = sqrtm(A)
  
  tn = function(e) {
    #norm(A_sqrt %*% e, "2")^2
    as.numeric(t(e) %*% A %*% e)
  }
  
  # tobs
  bR = R %*% bhat
  tobs = as.numeric(t(bR) %*% solve(R%*% V2b %*% t(R)) %*% (bR)) 
  # randomization distribution
  
  f0 = lm(y ~ z)
  e0 = f0$residuals
  n = length(e0)
  
  if(!use_restricted) {
    e0 = fit$residuals
  }
  
  tvals = replicate(RAND_REPS, {
    s = sign(rnorm(n))
    tn(s*sample(e0))  
    # tn(e*s)
  })
  
  p = mean(tvals >= tobs)
  (p <= 0.05)
}

# Main simulation
main_sim = function(nreps=100) {
  
  all_n = c(10, 25, 50, 100)
  # all_n = c(10, 20)
  cols = c("n_size","rri_lin", "rri_lin_R", "rri_quad", "rri_quad_R")
  Results = matrix(0, nrow=0, ncol=length(cols))
  colnames(Results) = cols
  
  summ_output = function(Rmat) {
    m = as.data.frame(Rmat)
    m1 = m %>% group_by(n_size) %>% summarize( rri_lin = 100*mean(rri_lin), 
                                               rri_lin_R = 100*mean(rri_lin_R), 
                                               rri_quad = 100*mean(rri_quad), 
                                               rri_quad_R = 100*mean(rri_quad_R),
                                               num=length(rri_lin))
    return(m1)
  }
  
  for(n in all_n) {
    for(j in 1:nreps) {
      x = runif(n, min=1, max=4)
      # x = rweibull(n, shape=0.3)
      z  = rnorm(n)
      eps = rnorm(n)*mean(abs(x))
      y = 0 + 0*x + z + eps
      
      l = RRI_lin(y, x, z, use_restricted = F)
      l_R = RRI_lin(y, x, z, use_restricted = T)
      q = RRI_quad(y, x, z, use_restricted = F)
      q_R = RRI_quad(y, x, z, use_restricted = T)
    
      Results = rbind(Results, c(n, l, l_R, q, q_R))
      m = summ_output(Results)
      print(m)
    }
    m = summ_output(Results)
    print(m)
  }
}

