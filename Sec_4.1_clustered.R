#' Panos Toulis (panos.toulis@chicagobooth.edu)
#' Invariant inference via residual randomization.
#' 
#' June, 2022
#' Reproduces Section 4.1 of paper, Table 1

rm(list=ls())
library(lmtest)
library(sandwich)
source("rri_functions.R")


#' Follows DGPs from Cameron 2008.
#' Y_{ig} = beta0 + beta1*X_ig + U_ig
#'    X_{ig} = Z_g + Z_{ig}
#'    U_ig = E_g + E_{ig}
#'    
DGP = function(J=10, Ng=30, het_err=F, Xg_scheme="norm", Eg_scheme="norm", Eig_scheme="norm", noise_ftr=1) {
  
  N = Ng*J # total units
  # h_errors = {T, F} to set heteroskedastic errors.
  #
  beta0 = 0
  beta1 = 0
  if(het_err) {
    beta0 = 1 # as in Cameron 2008.
    beta1 = 0
  }
  
  ## Covariates
  Zig =  matrix(rnorm(N), ncol=J)  # Zig, unit x cluster
  Zg = NA # cluster effects
  if(Xg_scheme=="norm") {
    Zg = matrix(rnorm(N), nrow=Ng, ncol=J, byrow=T)   
  } else if(Xg_scheme=="lognorm") {
    Zg = matrix(.46 * rlnorm(N), nrow=Ng, ncol=J, byrow=T)
  } else {
    stop(sprintf("%s X-scheme not implemented", Xg_scheme))
  }
  # Zg = cluster effect, Z_ig = idiosyncratic . Both unit x cluster. (Ng x J)
  Xig =  Zg + noise_ftr * Zig   # covariates
  
  # Errors
  Eg = NA
  if(Eg_scheme=="norm") {
    Eg = matrix(rnorm(N), nrow=Ng, ncol=J, byrow=T)
  } else if(Eg_scheme=="zero") {
    Eg = matrix(0, nrow=Ng, ncol=J, byrow=T)
  } else {
    stop(paste("Scheme", Eg_scheme, "not implemented"))
  }
  # Eg =  cluster error (Ng x J)
  
  # Vig = 1 if homoskedastic errors.
  Vig = het_err * (3 * abs(Xig)) + (1-het_err) * matrix(1, nrow= Ng, ncol=J)
  
  Eig = NA
  if(Eig_scheme=="norm") {
    Eig = Vig * matrix(rnorm(N), nrow=Ng) #
  } else if(Eig_scheme=="t") {
    Eig = Vig * matrix(rt(N, df=3), nrow=Ng) #
  }  else {
    stop(sprintf("%s Eig-scheme not implemented", Eig_scheme))
  }
 # errors
  Uig = Eg + Eig  # errors
  #
  x = as.vector(Xig)  # stack all clusters one-by-one., traverses x by column.
  # inverse: as.matrix(x, nrow=Ng)
  u = as.vector(Uig)
  
  y = beta0 + beta1 * x + u  # outcomes
  Yig = matrix(y, nrow=Ng)  # re-matrix a vector.
  
  ## we will save everything in a list
  DATA = list()
  DATA$x = x
  DATA$y = y
  DATA$Xg = list()
  DATA$Yg = list()
  # store variables per cluster.
  for(g in 1:J) {
    Xg = matrix(c(rep(1, Ng), Xig[, g]), ncol=2)
    Yg = matrix(Yig[, g], ncol=1)
    DATA$Xg[[g]] = Xg
    DATA$Yg[[g]] = Yg
  }
  DATA$clust = as.vector(matrix(1:J, nrow=Ng, ncol=J, byrow = T))
  DATA$J = J
  DATA$Ng = Ng
  
  return(DATA)
}

CHECK_DATA = function(DATA) {
  # Checks internal consistency of data.
  stopifnot(setequal(names(DATA), c("x", "Xg", "y", "Yg")))
  all(unlist(DATA$Yg) == DATA$y)
}

# OLS test and cluster-robust test.
#
ols_test = function(DATA, H0) {
  y = DATA$y; x = DATA$x
  fit = lm(y ~ x)
  ci = confint(fit)
  left = ci[2, 1]; right = ci[2, 2]
  (left > H0) || (right < H0)
}

# Cluster-robust test.
crve_test = function(DATA, H0) {
  f = lm(y ~ x, data=DATA)
  ci = coefci(f, vcov=vcovCL(f, cluster=DATA$clust))
  left = ci[2, 1]; right = ci[2, 2]
  (left > H0) || (right < H0)
}


single_sim = function(DATA, beta1_H0) {
  # Standard OLS and cluster-robust tests.
  ols_test = ols_test(DATA, H0 = beta1_H0)
  cr_test = crve_test(DATA, H0 = beta1_H0)
  
  # Residual randomization
  ## BUILD the clusters
  J = DATA$J
  cl0 = list()
  for(gg in 1:J) {
    cl0[[gg]] = which(DATA$clust==gg)
  }
  
  #' cl0 has the clustering we want
  #' Prepare the data for RR.
  y = DATA$y; X = cbind(1, DATA$x)
  rr.control = basic_rrcontrol
  rr.control$type = "sign"
  rr.control$clustering = cl0
  #
  sign_test = rrtest(y, X, a=c(0,1), a0=beta1_H0, rr.control)
  rr.control2 = rr.control
  rr.control2$type = "double"
  double_test = rrtest(y, X, a=c(0,1), a0=beta1_H0, rr.control2)
  
  return(data.frame(ols=ols_test, crve=cr_test, rr_sign=sign_test, rr_double=double_test))
}

#' Produces Table 1 for given J (number of clusters.)
main_sim = function(beta1_H0=0, J=10, nreps=5000) {
  
  methods = c("ols", "crve", "RR-sign", "RR-double")
  cols = c(methods, "het", "Xg", "Eg", "Eig", "H0")
  Results = matrix(0, nrow=0, ncol=length(cols))
  colnames(Results) = cols
  
  t0 = proc.time()[3] # start time.
  
  all_Xg_scheme = c("norm","lognorm")
  all_Eg_scheme = c("zero", "norm")
  all_Eig_scheme = c("norm", "t")
  all_HET = c(F, T)
  
  # Main iteration.
  for(het in all_HET) {
    for(i in 1:length(all_Xg_scheme)) {
      for(j in 1:length(all_Eg_scheme)) {
        for(k in 1:length(all_Eig_scheme)) {
          
          M = matrix(0, nrow=0, ncol=length(methods)) # interim results.
          colnames(M) = head(cols, length(methods))
          
          for(iter in 1:nreps) {
            # 1. generate data
            d = DGP(J=J, het_err = het,
                     Xg_scheme = all_Xg_scheme[i],
                     Eg_scheme = all_Eg_scheme[j],
                     Eig_scheme = all_Eig_scheme[k])
            
            test_results = single_sim(d, beta1_H0)
  
            # Store results.
            M = rbind(M, c(test_results$ols, test_results$crve, 
                           test_results$rr_sign, test_results$rr_double))
            t1 = proc.time()[3]
            if(t1 -t0 > runif(1, min=5, max=20)) {
              print(sprintf("iter = %d / %d --- H0 = %.2f, J=%d", iter, nreps, 0, J))
              print(round(100*colMeans(M), 3))
              t0 = t1
            }
          }
          
          # Simulation setting complete. Store to RESULTS
          cm = as.data.frame(t(round(100*colMeans(M), 3)))
          R1 = cbind(cm, data.frame(het=het, Xg=c(all_Xg_scheme[i]), Eg=c(all_Eg_scheme[j]), Eig=all_Eig_scheme[k]))
          Results = rbind(Results, R1)
          
          print("=--=-=-  Table 1 Results   ===-===-===-")
          print(Results)
          
          save(Results, file="Table1.rda")
        }
      }
    }
  }
}

#' Illustrates some key examples from Table 1.
#' Also presented in accompanying vignette.
quick_examples = function(nreps=999) {
  
  set.seed(210)
  methods = c("OLS", "ClusterOLS", "RR-sign", "RR-double")
  init_results = function() {
    R = matrix(0, nrow=0, ncol=length(methods))
    colnames(R) = methods
    return(R)
  }
  
  # Example 1: Failure of CR errors. Over-reject with log-normal X.
  Results = init_results()
  print("Example with log-normal Xg, no cluster effects.")
  for(i in 1:nreps) {
    d = DGP(J=10, het_err = FALSE, Xg_scheme = "lognorm", Eg_scheme = "zero", Eig_scheme = "norm")
    out = single_sim(d, beta1_H0=0)
    Results = rbind(Results, as.matrix(out))
    if(nrow(Results) %% 10 ==0) {
      print(paste("j=", nrow(Results), " / ", nreps))
      print(round(100*colMeans(Results), 3))
    }
  }
  
  # Example 2: Failure of CR errors wiht heteroskedasticity, log-normal X.
  Results = init_results()
  set.seed(233)
  print("Example with log-normal Xg, no cluster effects, heteroskedastic error.")
  for(i in 1:nreps) {
    d = DGP(J=10, het_err = TRUE, Xg_scheme = "lognorm", Eg_scheme = "zero", Eig_scheme = "norm")
    out = single_sim(d, beta1_H0=0)
    Results = rbind(Results, as.matrix(out))
    if(nrow(Results) %% 10 ==0) {
      print(paste("j=", nrow(Results), " / ", nreps))
      print(round(100*colMeans(Results), 3))
    }
  }
}