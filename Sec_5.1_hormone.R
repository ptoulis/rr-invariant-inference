#' Panos Toulis (panos.toulis@chicagobooth.edu)
#' Invariant inference via residual randomization.
#' 
#' June, 2022
#' Reproduces Section 5.1 of paper, Fig 1 and Table 4.
#' 
rm(list=ls())
## Load required Packages.
library(bootstrap)
library(ggplot2)
library(latex2exp)

source("rri_functions.R")

data("hormone")

# Bootstrap confidence interval
boot_ci = function(y, X, nBoot=2000) {
  i0 = 2 # coefficient to focus on (1=bias, 2=slope)
  bhat = coef(lm(y ~ X+0))[i0]
  print(bhat)
  n = length(y)
  se = sd(replicate(nBoot, {
    i = sample(1:n, replace=T)
    yi = y[i]; Xi = X[i, ]
    coef(lm(yi ~ Xi + 0))[i0]
  }))
  print(sprintf("Bootstrap se=%.4f", se))
  return(c(bhat - 2*se, bhat + 2*se))
}

y0 = hormone$amount
X0 = matrix(hormone$hrs, ncol=1); X0 = cbind(rep(1, nrow(X0)), X0)
print("Bootstrap interval")
print(boot_ci(y0, X0))


invert_rri_test = function(name, rrtest_fn)  {
  
  b1_vals = seq(-0.1, 0,  length.out=250)
  pvals = sapply(b1_vals, function(b0) {
    rrtest_fn(b0)
  })
  
  left = min(b1_vals[pvals >= 0.05])
  right = max(b1_vals[pvals >= 0.05])
  print(paste("95% CI for method", name))
  print(round(c(left, right), 4))
  
  return(list(b1_vals=b1_vals, pvals=pvals))
}


#' Produces Table 4.
main_sim = function() {
  # Main simulation.
  rr.control = basic_rrcontrol
  rr.control$return_pval = TRUE
  a = c(0, 1); 
  
  # perm test
  perm_test = function(beta1_H0) {
    rr.control$type = "perm"
    rrtest(y0, X0, a, a0=beta1_H0, rr.control=rr.control)
  }

    invert_rri_test("full perm", perm_test)
  
  # sign test
  sign_test = function(beta1_H0) {
    rr.control$type = "sign"
    rrtest(y0, X0, a, a0=beta1_H0, rr.control=rr.control)
  }

  invert_rri_test("full sign", sign_test)
  
  # perm cluster
  perm_clust_test = function(beta1_H0) {
    rr.control$type = "perm"
    rr.control$clustering = list(1:9, 10:18, 19:27) # by manufacturer
    rrtest(y0, X0, a, a0=beta1_H0, rr.control=rr.control)
  }
  
  invert_rri_test("cluster perm", perm_clust_test)
  
  # double cluster
  perm_clust_double = function(beta1_H0) {
    rr.control$type = "double"
    rr.control$clustering = list(1:9, 10:18, 19:27) # by manufacturer
    rrtest(y0, X0, a, a0=beta1_H0, rr.control=rr.control)
  }
  
  invert_rri_test("cluster double", perm_clust_double)
}

#' Produces Figure 1.
plot_pvals = function() {
  
  rr.control = basic_rrcontrol
  rr.control$return_pval = TRUE
  a = c(0, 1); 
  
  # perm test
  perm_test = function(beta1_H0) {
    rr.control$type = "perm"
    rrtest(y0, X0, a, a0=beta1_H0, rr.control=rr.control)
  }
  
  out = invert_rri_test("full perm", perm_test)
  pvals = out$pvals
  b1_vals = out$b1_vals
  
  i = which(pvals > 0.025); f = loess(pvals ~ b1_vals,span = .1)
  D = data.frame(beta1=b1_vals, pval=pvals)
  
  rLeft = min(b1_vals[i]); rRight = max(b1_vals[i])
  g = ggplot(data=D, aes(beta1, pval)); 
  g = g + geom_area(color="gray", fill="darkgray") + 
    stat_smooth(geom = 'area', method = 'loess', span = .1, alpha = .25, fill = "blue")
  g = g + geom_hline(yintercept = 0.025, lty="dashed", lwd=.8, col="red")
  g = g + geom_linerange(x=rLeft, ymin=0, ymax=.1, lwd=1.5, col="red")
  g = g + geom_linerange(x=rRight, ymin=0, ymax=.1, lwd=1.5, col="red")
  g = g + theme(text = element_text(size=25)) + xlab(TeX("$\\beta_1^0"))
  g = g + theme(axis.title.y = element_text(margin = margin(t = 0, r = 40, b = 0, l = 0)))
  g = g + theme(axis.title.x = element_text(margin = margin(t = 40, r = 0, b = 0, l = 0)))
  plot(g)
  
  # print("Randomization interval")
  print(sprintf("Randomiz. CI %.3f, %.3f", rLeft, rRight))
}
