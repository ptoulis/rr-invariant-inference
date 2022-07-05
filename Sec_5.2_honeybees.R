#' Panos Toulis (panos.toulis@chicagobooth.edu)
#' Invariant inference via residual randomization.
#' 
#' June, 2022
#' Reproduces Section 5.2 of paper. Honeybees example. Ch 19 in Zuuer at al textbook.
#' 
rm(list=ls())

library(nlme)
library(dbarts)
library(Matrix)

source("rri_functions.R")
load("Section_5.2_honeybees.rdata")

head(Bees)
# Author pre-processing
Bees$fHive = factor(Bees$hive)
Bees$Lspobee = log(Bees$spore_density + 1, base=10)
Bees$Infection01 = factor(Bees$infection > 0)


# authors model
summary(lme(Lspobee ~ Infection01, random=~1| fHive, 
            data=Bees, method="REML", 
            weights=varIdent(form=~1 | Infection01)))

# RRI - double invariance
# create clustering
hive_ids = unique(levels(Bees$fHive))
cl = list()
for(j in hive_ids) {
  cl[[length(cl)+1]] = which(Bees$fHive==j)
}
# cl = {(1,2,3), (4, 5, 6), }... LIST of groupings by hives.

rr.control = basic_rrcontrol
rr.control$type = "double"  # double invariance
rr.control$num_r = 5000
rr.control$clustering = cl # cluster by hive.
y = Bees$Lspobee; X = cbind(1, Bees$Infection01); a = c(0, 1); 
rr.control$return_pval = TRUE # need to invert the test.

# single H0 test
rtest = function(beta1_H0) {
  rrtest(y, X, a, beta1_H0, rr.control)
}

# Invert the test.
invert_rrtest(rtest, CI_initial = c(-5, 5))


