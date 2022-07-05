#' Panos Toulis (panos.toulis@chicagobooth.edu)
#' Invariant inference via residual randomization.
#' 
#' June, 2022
#' Reproduces Section 4.2 of paper, Table 2
rm(list=ls())

library(sandwich)
library(lmtest)
library(arm)
library(fixest)

source("rri_functions.R")

# Cluster robust variance estimation function
robust.se.nodfc = function(model, cluster){
  require(sandwich)
  require(lmtest)
  M <- length(unique(cluster))
  N <- length(cluster)
  K <- model$rank
  dfc <- 1
  uj <- apply(estfun(model), 2, function(x) tapply(x, cluster, sum))
  rcse.cov <- dfc * sandwich(model, meat = crossprod(uj)/N)
  rcse.se <- coeftest(model, rcse.cov)
  return(list(rcse.cov, rcse.se))
}

robust.se  = function(model, cluster){
  require(sandwich)
  require(lmtest)
  M <- length(unique(cluster))
  N <- length(cluster)
  K <- model$rank
  dfc <- (M/(M - 1)) * ((N - 1)/(N - K))
  uj <- apply(estfun(model), 2, function(x) tapply(x, cluster, sum))
  rcse.cov <- dfc * sandwich(model, meat = crossprod(uj)/N)
  rcse.se <- coeftest(model, rcse.cov)
  return(list(rcse.cov, rcse.se))
}


#################
# Coverage
#################
reject = function(l, r, beta1_H0) {
  (l > beta1_H0) || (r < beta1_H0)
}

## 
methods = c("HC2", "two-way", "random-effects", "dyadCL", "RR-dyadic")
cols = c(methods, "error_type", "x_type", "n")
Results = matrix(0, nrow=0, ncol=length(cols))
colnames(Results) = cols
error_types = c("norm", "lognorm")
x_types = c("norm", "lognorm")

na_counter = matrix(0, nrow=0, ncol=length(methods))

single_sim = function(dataUp, beta1_H0, only.methods=F) {
  # OLS
  if(only.methods) {
    return(c("HC2", "2way", "random-effects", "dyad-CL", "RR-dyadic"))
  }
  N = max(dataUp$dyad1)
  index = 1:N
  fit = lm(dY~dX, data=dataUp)
  bhat = coef(fit)
  
  # HC2
  ci = coefci(fit, vcov=vcovHC(fit, type="HC2"))
  left = ci[2, 1]; right = ci[2, 2]
  hc2 = reject(left, right, beta1_H0)
  
  # 2way-clustered
  dyad.mat = cbind(as.integer(dataUp$dyad1), as.integer(dataUp$dyad2))
  fit2 = lm(dY~dX + factor(dyad1) + factor(dyad2), data=dataUp)
  ci = coefci(fit2, vcov=vcovCL(fit2, cluster=~dyad1+dyad2, fix=F))
  left = ci[2, 1]; right = ci[2, 2]
  # dyad.fit = feols(dY~dX | dyad1 + dyad2, data = dataUp)
  #  s = summary(dyad.fit)
  #  ci = s$coefficients[1] + c(-1, 1)*1.96*s$se[1]
  #  left = ci[1]; right=ci[2]
  two_way = reject(left, right, beta1_H0)
  
  # legacy code: Does not work well.
  # Dyadic cluster robust via multiway decomposition
  # for(i in 1:length(index)){
  #   iUp <- index[i]
  #   clusUp <- apply(dyad.mat,
  #                   1,
  #                   function(x)as.numeric(iUp %in% x))
  #   clusIndexUp <- clusUp*-99 + (1-clusUp)*1:nrow(dyad.mat)
  #   if(i==1){dcrUp <- robust.se.nodfc(fit, clusIndexUp)[[1]]}
  #   if(i>1){dcrUp <- dcrUp + robust.se.nodfc(fit, clusIndexUp)[[1]]}
  # }
  # # substract naive CR:
  # dcrUp2 = dcrUp - robust.se.nodfc(fit, dataUp$dyads)[[1]]
  # Vhat = dcrUp2 - (length(index)-2)*vcovHC(fit,type="HC0")
  # sehat = sqrt(diag(Vhat))
  # left = bhat[2] - 1.96*sehat[2]
  # right = bhat[2] + 1.96*sehat[2]
  dyadCL = reject(left, right, beta1_H0)
  
  
  # Random effects
  # linear mixed model
  fit.re = lmer(dY~dX + (1|dyad1) + (1|dyad2), data=dataUp)
  bRE = fixef(fit.re)
  seRE = se.fixef(fit.re)
  
  left = bRE[2] - 1.96*seRE[2]
  right = bRE[2] + 1.96*seRE[2]
  reff = reject(left, right,  beta1_H0)
  
  # randomization
  y=dataUp$dY; X=cbind(1, dataUp$dX); a=c(0, 1); a0=beta1_H0
  rr.control = basic_rrcontrol
  rr.control$type = "dyadic"
  
  rr = rrtest(y,X,a,a0, rr.control)
  
  c(hc2, two_way, reff, dyadCL, rr)
}

DGP = function(N, error_type, x_type) {
  index = 1:N
  n = choose(N, 2)
  # Dyads
  dyads =  as.matrix(apply(t(combn(index, 2)), 1, function(x) paste(x[1], x[2],sep="-")))
  # dyads = "1-2", "1-3", ...."9-10"
  dyad.mat = t(apply(dyads, 1, function(x) as.numeric(unlist(strsplit(x,"-")))))
  # [,1] [,2]
  # [1,] "1"  "2" 
  # [2,] "1"  "3" 
  # [3,] "1"  "4" 
  #...
  # Simulation
  # Generate the data
  a = rnorm(N)  #  fixed effects.
  if(error_type=="lognorm") {
    a = 0.5*rlnorm(N)
  }
  X = rnorm(N)
  if(x_type=="lognorm") {
    X = 0.5*rlnorm(N)
  }
  ### 
  dX <- da <- NA
  for(i in 1:n) {
    da[i] =  sum(a[dyad.mat[i,]])  # a_i
    dX[i] = abs(diff(X[ dyad.mat[i,] ]))
  }
  
  # model
  # Y_ij = a_i + a_j + beta1*|X_i - X_j| + u_i
  # da_ij = e_i + e_j.
  dY = 0 + da + 1 * dX + rnorm(n)
  
  # Full data frame
  dataUp = data.frame(dyads, dY, dX, 
                      dyad1 = as.numeric(dyad.mat[,1]), 
                      dyad2 = as.numeric(dyad.mat[,2]))
  return(dataUp)
}

main_sim = function(beta1_H0=1, nreps=100) {
  
  for(N in c(20, 30)) {
    for(er_type in error_types) {
      for(x_type in x_types) {
        
        M = matrix(0, nrow=0, ncol=length(methods)) ## results.
        colnames(M) = methods
        t0 = proc.time()[3]
        
        for(irep in 1:nreps) {
          
          dataUp = DGP(N, er_type, x_type) # data
          results_irep = single_sim(dataUp, beta1_H0) # methods
          
          ## Store results.
          M = rbind(M,  results_irep)
          t1 = proc.time()[3]
          if(t1-t0 > 5) {
            print(paste("sim=",  irep, "/",  nreps, " N=", N, "n=", n))
            print(round(100*colMeans(M, na.rm=T),3))
            print("NA % ")
            print(round(100*colMeans(is.na(M)), 3))
            t0 = t1
          }
        }
        
        # STORE RESULTS
        m = round(100*colMeans(M, na.rm=T), 3)
        na_counter = rbind(na_counter, colMeans(is.na(M)))
        cm = as.data.frame(t(m))
        R1 = cbind(cm, data.frame(error=er_type, x_type=x_type, n=n))
        Results = rbind(Results, R1)
        print("-=-===--  RESULTS  -=-===--")
        print(Results)
        print("NAs")
        print(round(100*na_counter, 3))
        save(Results, file="Table2.rda")
      }
    }
  }
}

#' 
quick_examples = function(nreps=100) {
  set.seed(209)
  rej = 0 
  # Reproduces column (e, x)=(normal, lognormal), n=190
  methods = single_sim(d1, 1, TRUE)
  for(irep in 1:nreps) {
    d1 = DGP(20, "norm", "lognorm")
    rej = rej + single_sim(d1, 1)
    names(rej) = methods
    # print results.
    print(round(100*rej/irep, 3))
  }
}




