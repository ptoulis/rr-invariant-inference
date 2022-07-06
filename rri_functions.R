#' v1.1 05/15
#' v1.2 06/23
#' v1.3 06/25 -- added 2way invariance

#' Creates the invariance function f(e) -> ge
#' @param type Character in {"perm", "sign", "double"}
#' @param clustering LIST of sets of indices.
#' 
create_invariance = function(n, type, clustering=list(), ...) {
  stopifnot("Should supply n=data dimension"=is.numeric(n))
  stopifnot("Invariance type should be a character"=is.character(type))

  perm = function(u) {
    stopifnot(length(u)==n)
    sample(u)
  }
  symm = function(u) {
    stopifnot(length(u)==n)
    sign(rnorm(length(u)))*u
  }
  
  if(length(clustering)==0) {
    if(type=="perm") {
      return(perm)
    } else if(type=="sign") {
      return(symm)
    } else if(type=="double") {
      return(function(u) symm(perm(u)))
    } else if(type=="dyadic") {
      N = polyroot(c(-2*n, -1, 1))
      if(any(abs(Im(N)) > 1e-16)) stop("Could not solve n=N(N-1)/2.")
      N = Re(N)[which(Re(N) > 0)]
      stopifnot("Could not solve n=N(N-1)/2."=floor(N)==N)
      # n = N (N-1)/2
      index = 1:N
      # Dyads
      dyads =  as.matrix(apply(t(combn(index, 2)), 1, function(x) paste(x[1], x[2],sep="-")))
      # dyads = "1-2", "1-3", ...."9-10"
      dyad.mat = t(apply(dyads, 1, function(x) as.numeric(unlist(strsplit(x,"-")))))
      # all dyad pairs, no duplicates.
      to_m = function(v) {
        stopifnot(length(v)==N*(N-1)/2)
        V = matrix(0, N, N)
        for(i in 1:nrow(dyad.mat)) { 
          row = dyad.mat[i,];
          V[row[1], row[2]] = v[i]
        }
        V = V + t(V)
        return(V)
      }
      to_v = function(V) {
        stopifnot(all(dim(V) == c(N, N)))
        v1 = as.vector(t(upper.tri(V) * V))
        setdiff(v1, 0)
      }
      
      # CHECK
      r1 = runif(N*(N-1)/2)
      stopifnot("to_m/to_v are not inverse."=all( abs(to_v(to_m(r1)) - r1) < 1e-10 ))
      
      fn = function(u) {
        stopifnot(length(u)==n)
        U = to_m(u) # into NxN matrix
        i = sample(1:N)
        Unew = U[i, i]
        u_new = to_v(Unew)
        return(u_new)
      }
    } else if(type=="2way") {
      
      arguments = list(...) 
      stopifnot("Should supply row-column dimension"=all(c("r_dim", "c_dim") %in% names(arguments)))
      
      r_dim = arguments$r_dim
      c_dim = arguments$c_dim
      stopifnot("R x C should be equal to n"=((r_dim*c_dim)==n))
      
      fn = function(u) {
        stopifnot(length(u)==n)
        U = matrix(u, nrow=dim_r, ncol=dim_c) # into Nr x Nc matrix
        
        Unew = U[sample(1:dim_r), sample(1:dim_c)]
        u_new = as.vector(Unew) # convention
        return(u_new)
      }
    }
    else {
      stop(paste("Not supported invariance = ", type))
    }
  } else {
    # Pre-process clustering
    # create Cluster indicator
    n1 = sum(unlist(lapply(clustering, function(A) length(A)))) # total units.
    stopifnot("Size of clustering should match n"=n1==n)
    C = rep(0, n1)
    for(j in 1:length(clustering)) {
      I_j = clustering[[j]]
      C[I_j] = j
    }
    # C = (1, 1, 1, 2, 2, 3, 3, 3, 3, ....) cluster indicator per unit.
    clust_perm = function(u) {
      stopifnot("Size of input should match size of clustering"=length(u)==n1)
      u1 = u
      for(j in 1:length(clustering)) {
        I_j = clustering[[j]]
        if(length(I_j) > 1) {
          u1[I_j] = sample(u[I_j])
        }
      }
      return(u1)
    }
   
    clust_symm = function(u) {
      stopifnot("Size of input should match size of clustering"=length(u)==n1)
      s_c = sign(rnorm(length(clustering))) # sample cluster signs
      s_all = s_c[C]
      return(s_all*u)
    }
    # proper structure of clustering
    if(type=="perm") {
      return(clust_perm)
    } else if(type=="sign") {
      return(clust_symm)
    } else if(type=="double") {
      return(function(u) clust_symm(clust_perm(u)))
    } else {
      stop(paste("Not supported invariance = ", type))
    }
  }
}

randomization_pval = function(tobs, tvals, two.sided=FALSE) {
  # tvals are additional draws
  N = length(tvals)
  p1 = (1/(N+1)) *(sum(tvals > tobs) + runif(1)*(1+sum(abs(tvals - tobs) < 1e-15)))
  if(!two.sided) {
    return(p1)
  } else {
    p2 = 2*min(p1, 1-p1)
    return(p2)
  }
}

critical_value = function(tvals, level) {
  u = unique(sort(tvals))
  v = sapply(u, function(s) sum(tvals <= s)) / length(tvals)
  # v = quantiles for each of order statistic of t_vals.
  u[min(which(v >= 1-level))]
}

#' The critical value implementation is conservative whenever t_vals has small support.
#' #' Let M = #values, and x_ = #values<c_a, x0 = #values=c_a and x+ = #values > c_a 
#' Then,  P(reject) <= a and P(reject) > a-x0/M. So, it depends on the multiplicity at M.
#' e.g., when M=30,  P(reject) > 0.05-1/30 = 0.01667. 
#' 
#' Using the randomization_pval() randomizes the decision to be close to the nominal level.
one_sided_test = function(tobs, tvals, level=0.05, return_pval=FALSE) {
  # c_a = critical_value(tvals, level=level)
  # phi = (tobs > c_a) # decision
  # phi     
  p1 = randomization_pval(tobs, tvals, two.sided=FALSE)
  if(return_pval) {
    return(p1)
  } else {
    return(p1 <= level)
  }
}

two_sided_test = function(tobs, tvals, level=0.05, return_pval=FALSE) {
  # c_1 = critical_value(tvals, level=level/2)
  # c_2 = critical_value(-tvals, level=level/2)
  # 
  # phi = (tobs > c_1) + (-tobs > c_2) # decision
  # phi > 0
  p2 = randomization_pval(tobs, tvals, two.sided=TRUE)
  if(return_pval) {
    return(p2)
  } else {
    return(p2 <= level)
  }
  
}

#' Transforms a vector of discrete values to clustering.
#' e.g., (A, A, B, B, A, C, C, C) -> LIST [[ {1,2,5}, {3,4}, {6,7,8} ]]
#' @param v Vector of discrete values.
vector_to_clustering = function(v) {
  u = unique(v)
  cl = list()
  for(u_class in u) {
    cl[[length(cl)+1]] = which(v==u_class)
  }
  return(cl)
}

basic_rrcontrol = list(type="perm",
                       clustering=list(),
                       g_invar=NULL,
                       test_type="two-sided",
                       level=0.05,
                       num_r=5000, 
                       return_pval=FALSE) 
#' Tests linear hypothesis
#' @param y Response
#' @param X nxp covariates. First column should be 1.
#' @param a p vector. 
#' @param a0 Scalar.We are testing H0: a'beta = a0. 
#'           We test beta_0 + beta_1 + ..beta_p =0 by default
#' @param control LIST that controls test.
rrtest = function(y, X,  a=rep(1, ncol(X)), a0=0, rr.control=basic_rrcontrol) {
  
  stopifnot("Data dimensions in y,X should match"=length(y)==nrow(X))
  n = length(y) # data dimension
  if(!all(X[,1]==1)) stop("No intercept column.")
  stopifnot("num of datapoints conflict."=nrow(X)==n)
  fit0 = .lm.fit(X, y)
  b_ols = coef(fit0)  # OLS estimator
  
  stopifnot(length(a)==ncol(X))
  # observed statistic
  Tn = sum(a*b_ols) - a0 # a'b^-a0
  
  # TODO: Speed up?
  S = solve(t(X) %*% X)
  aSa = sum(t(a) %*% S %*% a)
  aQ = as.vector(t(a) %*% S %*% t(X))
  stopifnot(length(aQ)==n)
  
  # restricted OLS -- TODO:Speed up?
  b_rls = b_ols - (Tn/aSa)* as.numeric(S %*% a)
  e_rls = as.numeric(y - X %*% b_rls) # restricted residuals
  
  # tn() function
  tn = function(e) { sum(aQ*e) }
  
  # Create invariance function
  g = NA
  if(is.null(rr.control$g_invar)) {
    g = create_invariance(n, rr.control$type, rr.control$clustering)
  } else {
    # print("[INFO] Using g() invariance function.")
    g = rr.control$g_invar
  }
  
  # Residual randomization test. Additional draws.
  tvals = replicate(rr.control$num_r, {
    tn( g(e_rls) )
  })
  # tvals = c(Tn, tvals)
  
  if(rr.control$test_type=="one-sided") {
    return(one_sided_test(tobs=Tn, tvals=tvals, level=rr.control$level, return_pval = rr.control$return_pval))
  } else {
    return(two_sided_test(tobs=Tn, tvals=tvals, level=rr.control$level, return_pval = rr.control$return_pval))
  }
}

#' Inverts a randomization test.
#' @param rrtest_fn Function that tests single null. MUST return a p-value.
#' 
invert_rrtest = function(rrtest_fn, CI_initial=c(-50,50), level=0.05, verbose=TRUE, tol=1e-4) {
  
  CI = CI_initial
  lower_prev = CI[1]
  upper_prev = CI[2]
  mid_prev = NA
  KEEP=T
  
  t0 = proc.time()[3]
  print(paste("[INFO]: Inverting randomization test.."))
  while(KEEP) {
    mid  = mean(CI) # (left+right)/2
    mid_prev = (lower_prev + upper_prev)/2
    p = rrtest_fn(CI[1]) # p-value
    if(p < level) {
      # We just rejected H0: lowerCI
      lower_prev = CI[1]
      CI[1] = (CI[1] + mid)/2 # increase lower endpoint
    } else {
      CI[1] = (lower_prev + CI[1])/2 # decrease lower
    }
  
    p = rrtest_fn(CI[2])  # p-value
    if(p < level) {
      # We just rejected H0: upperCI
      upper_prev = CI[2]
      CI[2] = (CI[2] + mid)/2 # decrease higher endpoint
    } else {
      CI[2] = (upper_prev + CI[2])/2 # increase higher endpoint
    }
    
    if(abs(lower_prev - CI[1]) < tol) { KEEP = FALSE }
    t1 = proc.time()[3]
    if(t1 - t0 > 5) {
      print(paste("-[INFO]: Current ", round(100*(1-level), 1), "% CI = [", round(CI[1], 3), ",", round(CI[2], 3), "]"))
      t0 = t1
    }
  }
  return(CI)
}

#' CI for particular beta component.
rrCI_param = function(y, X, param, rr.control, level=0.05) {
  stopifnot(param > 0 & param <= ncol(X))
  a_j = rep(0, ncol(X)); a_j[param] = 1
  rtest = function(b) {
    rrtest(y, X,a=a_j, a0=b, rr.control=rr.control)
  }
  invert_rrtest(rtest, level=level)
}