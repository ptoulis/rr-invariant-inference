#' Panos Toulis (panos.toulis@chicagobooth.edu)
#' Invariant inference via residual randomization.
#' 
#' June, 2022
#' Reproduces Section 4.3 of paper, Table 3

rm(list=ls())

source("rri_functions.R")
require(sandwich)
require(plyr)
require(mvtnorm)

n = 30 # total units
N1 = 3; N0 = n-N1  # treated and control units.

#' Calculates all the cluster sign transformations.
#' Code is old. Can be improved. Look at create_invariance() in rri_functions().
inv_Group = function(ds) {
  y = ds$y
  D = ds$d
  stopifnot(mean(D)==.1);
  stopifnot(length(y)==length(D))
  ## produce exact test
  
  i1 = which(D==1); i0 = which(D==0)
  N1 = length(i1); N0 = n-N1
  
  # Build the solutions.
  mgen = function(a) matrix(c(a, a, -a, a, -a, a, -a, a, a), nrow=3, byrow=T)
  
  M1 = mgen(1)
  gen_row = function(a, j) {
    stopifnot(j<7)
    if(j<4) return(mgen(a)[j,])
    -mgen(a)[j-3,]
  }
    
  M0 = matrix(0, nrow=0, ncol=n)
  for(i in 1:6) {
    s1 = gen_row(1, i)
    jrange = c(i) # 1:3 # 3 * (i>3) + 1:3
    for(j in jrange) {
      s0 = gen_row(rep(1, 9), j)
      S = rep(0, n); S[i1] = s1; S[i0] = s0
      M0 = rbind(M0, S)
    }
  }
  M0 = rbind(M0, rep(1, n))
  M0 = rbind(M0, rep(-1, n))
  rownames(M0) = NULL
  
  ## M0 is specil matrix.
  # Check treated units.
  stopifnot(setequal(rowSums(M0[, i1]), c(-3, -1, 1, 3)))
  # Check control units.
  stopifnot(setequal(rowSums(M0[, i0]), c(-27, -9, 9, 27)))
  # Check closure.
  for(i in 1:nrow(M0)) {
    for(j in 1:nrow(M0)) {
      x1 = M0[i, ]; x2 = M0[j, ]
      a = apply(M0, 1, function(row) sum(abs(row - x1 * x2)))
      # print(a)
      stopifnot(sum(a==0) > 0)
    }
  }
  rownames(M0) = NULL
  return(M0)
}

#' For use with rri_functions() if necessary
create_C0_clustering = function(ds) {
  d = ds$d
  stopifnot(all(d %in% c(0, 1)))
  i1 = which(d==1)
  i0 = which(d==0)
  stopifnot("Control/Treated should be integer"=(length(i0) %%  length(i1) == 0))
  K = as.integer(length(i0) / length(i1)) # controls / treated unit
  cl = list() 
  while(length(i1) > 0) {
    treated = head(i1, 1)
    controls = sample(i0, size=K, replace=F)
    cl[[length(cl)+1]] = c(treated, controls)
    i1 = setdiff(i1, treated)
    i0 = setdiff(i0, controls)
  }
  return(cl)
}

# Exact test.
test_exact = function(ds, beta1_H0) {
  #
  fit = lm(y ~ d, data=ds)
  bhat = coef(fit)
  Tn = bhat[2] - beta1_H0
  # restricted residuals
  e_r = ds$y - mean(ds$y - beta1_H0*ds$d) - beta1_H0*ds$d
  
  X = cbind(1, ds$d)
  tn = function(u) coef(lm(u ~ X + 0))[2]  # tn(eps) = Tn if H0 is true
  
  # Test begins
  M = inv_Group(ds)   # tn(M_i' eps) = tn(M_i' e_r) => exact test.
  I1 = which(rowSums(M)==n)
  tvals = sapply(setdiff(1:nrow(M), I1), function(i) tn(M[i,]*e_r))
  two_sided_test(Tn, tvals)
}

err_mixture = function(nsampl) {
  # a = sample(c(-1, 1), size=nsampl, T)* sample(log(1:nsampl))
  # a/sd(a)
  B  = rbinom(nsampl, size=1, prob=0.5)
  a = B * rnorm(nsampl, mean=-5, sd=0.5) + (1-B) * rnorm(nsampl, 5, 0.5)
  a/5
}

# DGP in this simulation
DGP = function(error_type="norm", sigma0=1) {
  D = sample(c(rep(1, N1), rep(0, N0)))  ## only 3 treated.
  #
  i1 = which(D==1); i0 = which(D==0)
  # N1 = length(i1); N0 = length(i0)
  # errors.
  eps1 = eps0 = NA
  
  if(error_type=="norm") {
    eps1 = rnorm(n)
    eps0 = sigma0 *rnorm(n)
  } else if(error_type=="t") {
    eps1 = rt(n, df=3)
    eps0 = sigma0 * rt(n, df=3)
  } else if(error_type=="cauchy") {
    eps1 = rcauchy(n)
    eps0 = sigma0*rcauchy(n)
  } else if(error_type=="mixture") {
    eps1 = err_mixture(n)
    eps0 = sigma0*err_mixture(n)
  } else {
    stop("Not supported error type.")
  }
  eps = D * eps1 + (1-D) * eps0
  # Sample data.
  y = -1  + 1 * D + eps
  return(data.frame(y=y, d=D))
}

#################
# Coverage
#################
reject = function(l, r, beta1_H0) {
  (l > beta1_H0) || (r < beta1_H0)
}

#' One single simulation over dataset "ds"
#' 
single_sim = function(beta1_H0, ds) {
  # 1. ols, BM correction
  fit = lm(y ~ d, data=ds)
  est = coef(fit)[2]
  se = sqrt(vcovHC(fit, type = "HC2")[2, 2])
  K_BM = (N0 + N1)^2 * (N0-1)*(N1-1) / (N1^2* (N1-1) + N0^2*(N0-1))
  se_BM = se * (qt(.975, df=K_BM)/1.96)
  left = est - 2*se_BM; right = est + 2*se_BM
  rej_BM = reject(left, right, beta1_H0)
  
  # 2. Randomization tests.
  # randomization
  y = ds$y; X = cbind(1, ds$d); a=c(0, 1); a0=beta1_H0
  rr.control = basic_rrcontrol
  rr.control$num.r = 5000
  rr.control$type = "sign"
  rr.control$clustering = as.list(1:n)
  
  # wild bootstrap (ish)
  wild = rrtest(y, X, a, a0, rr.control=rr.control)
  
  # Exact randomization
  # rr.control$clustering = create_C0_clustering(ds)
  # exact = rrtest(y, X, a, a0, rr.control=rr.control)
  exact = test_exact(ds, beta1_H0)

  v = c(rej_BM, wild, exact)
  names(v) = c("BM", "wild boot", "exact")
  return(v)
}
# Main simulation for Behrens-Fisher problem.
#
# We consider (1) normal errors; (2) symmetric, but not normal errors.
#
main_sim = function(beta1_H0, sigma0=c(.5, 1, 2, 5, 10), nreps=100) {
  ## 
  methods = c("BM", "RR-sign", "RR-exact")
  cols = c(methods, "H0", "sigma_0", "tr_pct" ,"error_type")
  Results = matrix(0, nrow=0, ncol=length(cols))
  colnames(Results) = cols
  err_type = c("norm", "t", "cauchy", "mixture")
  
  for(ii in 1:length(err_type)) {
    for(sig0 in sigma0) {
      ## Iteration paramers: %treated, error type, sigma0
      t0 = proc.time()[3]
      M = matrix(0, nrow=0, ncol=length(methods))  # Only keep the numeric values here.
      colnames(M) = methods
      ## Main iteration
      for(iter in 1:nreps) {
        # Covariates.
        
        dat = DGP(err_type[ii], sig0) # data 
        iter_results = single_sim(beta1_H0, dat)
        
        # Store results.
        M = rbind(M, iter_results) # , sig0, ii, beta1_H0, round(100*ntreat/n, 2))) 
        
        t1 = proc.time()[3]
        if(t1 - t0 > runif(1, min=5, max=10)) {
          print(round(100*colMeans(M), 3))
          print(sprintf("Error type=%s, sigma0=%.2f, H0=%.2f, true=%.2f, N1=%d, n=%d, i=%d/%d -- ", 
                        err_type[ii], sig0, beta1_H0, 1, N1, n, iter, nreps))
          t0 = t1
        }
      }
      ## iterations done
      cm = as.data.frame(t(round(100*colMeans(M), 3)))
      R1 = cbind(cm, data.frame(H0=beta1_H0, sigma_0=sig0, tr_pct=100*N1/n, error_type=err_type[ii]))
      Results = rbind(Results, R1)
      #
      print("--===--     RESULTS    - --=-=-===-")
      print(Results)
    }
  }
}


## Run simulation
quick_examples  = function() {
  rej = 0
  set.seed(237)
  for(j in 1:1000) {
    ds = DGP("mixture", sigma0 = 3)
    v = single_sim(beta1_H0 = 1, ds)
    rej = rej + v
    if(j %% 100==0) {
      print(paste("iter=", j))
      print(round(100*(rej/j), 3))
    }
  }
}
