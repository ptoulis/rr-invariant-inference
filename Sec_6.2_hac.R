#' Panos Toulis (panos.toulis@chicagobooth.edu)
#' Invariant inference via residual randomization.
#' 
#' June, 2022
#' Reproduces Section 6.2 of paper, Table 7
#' rm(list=ls())
library(sandwich)
library(lmtest)
source("rri_functions.R")

rho = 0.8
# Other params
n = 100 # sample size
all = expand.grid(1:n, 1:n)
A_rho = apply(expand.grid(1:n, 1:n), 1, function(row) rho^abs(row[1] - row[2]))
A_rho = matrix(A_rho, nrow=n)
A_rho = lower.tri(A_rho, diag = T) * A_rho

err_mixture = function(nsampl) {
  # a = sample(c(-1, 1), size=nsampl, T)* sample(log(1:nsampl))
  # a/sd(a)
  B = rbinom(nsampl, size=1, prob=0.5)
  a = B * rnorm(nsampl, mean=-5, sd=0.5) + (1-B) * rnorm(nsampl, 5, 0.5)
  a/5
}

DGP = function(x_scheme="norm", e_scheme="norm") {
  
  u = x = NA
  if(e_scheme=="norm") {
    eps = rnorm(n)
  } else if (e_scheme=="ac-norm") {
    eps = A_rho %*% rnorm(n)
  } else if (e_scheme=="mixture") {
    eps = err_mixture(n)
  } else if (e_scheme=="ac-mixture") {
    eps = A_rho %*% err_mixture(n)
  } else {
    stop(e_scheme)
  }
  
  if(x_scheme=="norm") {
    x = rnorm(n)
  } else if (x_scheme=="ac-norm") {
    x = A_rho %*% rnorm(n)
  } else if (x_scheme=="lognorm") {
    x = rlnorm(n)
  } else if (x_scheme=="ac-lognorm") {
    x = A_rho %*% rlnorm(n)
  } else {
    stop(x_scheme)
  }
  
  y = 1 + 0 * x + eps
  return(data.frame(Yt=y, Xt=x))
}

create_Cn_plus = function(e) {
  z = abs(diff(sign(e)))
  cl = list()
  cl[[1]] = 1
  for(j in 2:length(e)) {
    if(z[j-1]==0) {
      # same sign as cluster
      cl[[length(cl)]] = c(cl[[length(cl)]], j)
    } else {
      # new cluster
      cl[[length(cl) +1]] = c(j)
    }
  }
  return(cl)
}

CHECK_Cn_plus = function() {
  e = rep(0, 100)
  for(i in 2:100) {
    e[i] = 0.8*e[i-1] + rnorm(1)
  }
  cl = create_Cn_plus(e)
  
  sapply(1:length(cl), function(i) length(cl[[i]]))

  stopifnot(all(sapply(1:length(cl), function(i) sd(sign(e[cl[[i]]]))==0), na.rm=T))
}

reject = function(left, right, beta1_H0) {
  (left > beta1_H0) || (beta1_H0 > right)
}

#' Single simulation
single_sim = function(beta1_H0, ds) {

  y = ds$Yt
  x = ds$Xt
  # OLS
  fit = lm(y ~ x)
  b0 = mean(y - beta1_H0 * x)
  eps_r = y - beta1_H0 * x - b0 # restricted residuals
  
  ## OLS test.
  ci = confint(fit)[2,]
  left = ci[1]
  right = ci[2]
  ols = reject(left, right, beta1_H0)
  
  ## HAC test
  ci = coefci(fit, df=Inf, vcov=vcovHAC)[2,]
  left = ci[1]
  right = ci[2]
  hac = reject(left, right, beta1_H0)
  
  # Randomization test.
  cl = create_Cn_plus(eps_r)
  K = max(sapply(1:length(cl), function(i) length(cl[[i]])))
  
  randz = 0; cond_randz = NA
  y = ds$Yt; # de-mean the y vector
  X = cbind(rep(1, n), ds$Xt)
  rand_test_run = 0
  if(K >= 3) {
    rr.control = basic_rrcontrol
    rr.control$type = "sign"
    rr.control$clustering = create_Cn_plus(eps_r)
    randz = rrtest(y,X, a=c(0,1), a0=beta1_H0, rr.control=rr.control)
    
    cond_randz = randz
    rand_test_run = 1
  }
  
  ret = c(ols, hac, randz, cond_randz, rand_test_run) # last two don't matter
  names(ret) = c("OLS", "HAC", "RR-reflection", "aux1", "aux2")
  return(ret)
}

#' Main simulation
main_sim = function(beta1_H0, nreps=100) {
  cols = c("ols", "hac",  "randomiz.", "randz_Kcond", "util(0,1)", "x_error", "e_error")
  Results = matrix(0, nrow=0, ncol=length(cols))
  colnames(Results) = cols
  t0 = proc.time()[3]
  
  all_x =  c("norm" , "lognorm", "ac-norm", "ac-lognorm")
  all_e = c("ac-norm", "ac-mixture")
  
  for(esc in all_e) {
    for(xsc in all_x) {
      
      numCols = cols[1:5]
      M = matrix(0, nrow=0, ncol=length(numCols))
      colnames(M) = numCols
      rcounter = 0
      
      for(i in 1:nreps) {
        
        ds = DGP(x_scheme = xsc, e_scheme = esc) # data
        results_i = single_sim(beta1_H0, ds) # methods/results
        rcounter = rcounter + tail(results_i, 1)
        
        M = rbind(M, results_i) # store results
        
        t1 = proc.time()[3]
        if(t1 - t0 > runif(1, min=5, max=10)) {
          print(sprintf("(%d/%d) -- H0=%.2f, beta0=%.2f, beta1=%.2f, Randomiz. util %% = %.2f%% (=%d), rho=%.2f", 
                        i, nreps, beta1_H0, 1, 0, 100 * rcounter/i, rcounter, rho))
          print(round(100*colMeans(M, na.rm = T), 3))
          t0  = t1
        }
      }
      
      cm = as.data.frame(t(round(100*colMeans(M, na.rm=T), 3)))
      R1 = cbind(cm, data.frame(x_error=xsc, e_error=esc))
      Results = rbind(Results, R1)
      #
      print("--===--     RESULTS    - --=-=-===-")
      print(Results)
    }
  }
  return(Results)
}

quick_examples = function() {
  set.seed(182)
  # Set rho=0.8
  R = 0
  for(irep in 1:1000) {
    ds = DGP(x_scheme = "ac-norm",e_scheme = "ac-norm")
    R = R + single_sim(beta1_H0=0, ds)
    if(irep %% 50==0) {
      print(paste("samples=",irep, "Rejection rates %="))
      print(round(100*(R[1:3]/irep), 3))
    }
  }
}
