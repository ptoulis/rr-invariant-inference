# Invariant Inference via Residual Randomization

Consider a simple linear model:
```
y = β0 + x * β1 + ε
```
where y=outcomes, Χ=covariates, ε=errors, and β=parameters.
In invariant inference, we only make structural invariance assumptions on the errors, but the data may be non-iid.
This type of inference is called *invariant inference*. 

Overall, invariant inference under residual randomization seems to be more robust in finite samples than robust OLS errors. Moreover, it is easier to apply in practice: all is required is to "plug-and-play" an appropriate invariance structure in the randomization procedure. See the paper for more details.

## Example 1. Simple exchangeable errors.
Let's illustrate with an example. First, we sample the data.
```
rm(list=ls())
source("rri_functions.R") # loads main library

set.seed(123)
n = 200
x = runif(n)
y = -1 + 0*x + rnorm(n)
```
We want to test H0: β1=0. This is a true null.
Before we run the residual randomization (RR), we first initialize its controls:
```
rr = basic_rrcontrol
rr$return_pval = TRUE
rr$g_invar = create_invariance(n, type="perm")
X = cbind(1, x)
```
The `rr` object specifies that we want to do inference under exchangeable errors (as discussed above).
We are now ready to run the basic RR test.
```
rrtest(y, X, a=c(0,1), a0=0, rr.control=rr) # tests a'beta = a0, pval = 0.59
```
The p-value is much larger than 0.05. We cannot reject the null.
Now, to check the power, let's test H0: β1=1 (false).
```
rrtest(y, X, a=c(0,1), a0=1, rr.control=rr) # pval = 2e-4  -> reject at 5%
 
```
We see that RR correctly rejects this null.
Finally, we can invert the test and produce confidence intervals:
```
rrCI_param(y, X, param=2, rr) # -0.6220583  0.3387069
confint(lm(y ~ x)) #   -0.6318838  0.3522660
```
We see that the RR-based confidence interval is very similar to OLS.

## Example 2. Clustered errors
We tweak the above code to produce clustered data.
```
rm(list=ls())
source("rri_functions.R")
set.seed(123)
n = 200
x = sort(runif(n))
I = c(rep(0, n*0.9), rep(1, 0.1*n)) # cluster indicator
sigma = (1-I)*0.1+ I*5 # larger variance when I=0
y = -1 + 0*x + rnorm(n, sd=sigma)
```
Use `plot(x, y)` to see the effects of this code. Roughly speaking, when x>0.9 the variance of y is much higher. This means that the data are not iid, but are heteroskedastic. Of course, standard OLS is off:
```
confint(lm(y ~ x)) #  0.4317599  2.0626657
```
OLS says the slope coefficient is strongly positive.
How can we account for such heteroskedasticity? It is easy with RR by considering the appropriate error structure. 
Here, the errors are exchangeable within the clusters specified by `I`. We then use the code:
```
rr$type = "perm"   # not necessary here because "perm" is set by default.
rr$clustering = vector_to_clustering(I)
rrtest(y,X, a=c(0, 1), a0=0, rr.control=rr)  # p-value = 0.09282712
```
We see that the RR procedure with the correct structure does not reject the null. 
Notice that we only changed the structure in the input of `rrtest`. The procedure remains exactly the same.

## Example 3. Wait, can't I just use robust errors?
Yes, you can, but these procedures tend to underperform in many practical settings, especially with small samples or heavy-tailed data.
Let's run a head-to-head comparison of RR against a "heteroskedasticity robust" method for illustration.
First, we write our DGP:
```rm(list=ls())
source("rri_functions.R")
library(lmtest)
set.seed(123)

DGP = function(x.uniform=TRUE) {
  n = 100
  x = 10*runif(n)
  if(!x.uniform) {
    x = 10*rweibull(n, shape=0.5)
  }
  y = -1 + 0*x + sqrt(abs(x))*rnorm(n) 
  data.frame(y=y, x=x)
}
```
This includes a choice to sample Weibull covariates, which are heavy-tailed.
A robust OLS test can be defined as follows:
```
ols_test = function(d) {
  fit = lm(y ~ x, data=d)
  V = vcovHC(fit, type = "HC2")
  ci = coefci(fit, df = Inf, vcov = vcovHC(fit, type = "HC0"))
  ci[2,1] > 0 | ci[2,2] < 0  # tests β1=0
}
```
This is using "HC2" robust standard errors.
The analogous RR test can be defined as follows:
```
rr = basic_rrcontrol
rr$type= "sign"  # we assume sign symmetric errors. Useful under general heteroskedasticity.
RR2_test = function(d) {
  X = cbind(1, d$x)
  y = d$y
  rrtest(y, X, a=c(0, 1), a0=0, rr.control=rr) # tests β1=0
}
```
Now, we are ready to kick off our simulation:
```
nreps = 500
rej = t(replicate(nreps, {
  d = DGP(x.uniform = TRUE)
  c(ols_test(d), RR2_test(d))
}))
colnames(rej) = c("OLS", "RR")
print(round(100*colMeans(rej), 3))

[1] OLS  RR 
    4.8 4.6 
```
In that first simulation, everything works ok for both methods. Our x are normal, our errors are normal, and we have 200 samples to work with. 
This is "simulationland" (favored in academic papers) where OLS with "robust errors" has no problems. 

In the next piece of code, we use heavy-tailed covariates, and the results change drastically:
```
rej = t(replicate(nreps, {
  d = DGP(x.uniform = FALSE)
  c(ols_test(d), RR2_test(d))
}))
colnames(rej) = c("OLS", "RR")
print(round(100*colMeans(rej), 3))

[1] OLS  RR 
    20   6 
```
These "robust OLS errors" are robust in name only! The randomization-based method here performs graciously, over-rejecting only slightly compared to OLS.


