#' Panos Toulis (panos.toulis@chicagobooth.edu)
#' Invariant inference via residual randomization.
#' 
#' June, 2022
#' Reproduces Section 5.3 of paper, Table 5.
#' 
rm(list=ls())
par(mfrow=c(1, 1))
library(haven)
library(tidyverse)
library(sandwich)
library(countrycode)

source("rri_functions.R")

get_continent = function(ccode) {
  ccode = as.character(ccode)
  a = subset(codelist, ioc==ccode | genc3c==ccode | genc2c==ccode | iso2c==ccode)
  as.character(a$continent)
}
get_countries = function(df) {
  as.character(unique(c(df$ccode1, df$ccode2)))
}

data_fname = "Sec_5.3_RoseEngel.RData"
#' Code to preprocess data. Already pre-processed in "_RData" file
#'
# 
load_data = function(remove_dup=FALSE) {
  print("Loading data...")
  trade = read_dta("rose_engel_trade/dat9512.dta")
  trade = trade %>% mutate(lrgdp = log((gdppc1*gdppc2)*(pop1*pop2)),
                           lrgdpcc = log((gdppc1*gdppc2)))

  # Reproduces exactly the regressions in pg1.log. "gravity model")
  summary(lm(lvalue ~ cu + ldist + lrgdpcc + lrgdp + regional + comlang + border, data=trade,
             na.action = na.omit))

  # Remove NA explicitly
  df = subset(trade, select=c(ccode1, ccode2, lvalue, cu, ldist, lrgdpcc, lrgdp, regional, comlang, border))
                             # gdppc1, gdppc2))
  df = na.omit(df)  # 4618 data pairs, 126 countries
  fit = lm(lvalue ~ cu + ldist + lrgdpcc + lrgdp + regional + comlang + border, data=df)
  summary(fit)

  if(remove_dup) {
    print("Removing duplicates..")
    # remove duplicates
    df = df %>% group_by(ccode1, ccode2) %>% summarize(ccode1=ccode1[1], ccode2=ccode2[1],
                                                       lvalue=mean(lvalue), cu=mean(cu), ldist=mean(ldist),
                                                       lrgdpcc=mean(lrgdpcc), lrgdp=mean(lrgdp),
                                                       regional=mean(regional), comlang=mean(comlang), border=mean(border))

  }
  # add continent
  print("Adding continent...")
  a  = sapply(df$ccode1, get_continent)
  df$cont1 = as.factor(as.character(a))
  b = sapply(df$ccode2, get_continent)
  df$cont2 = as.factor(as.character(b))

  return(df)
}


## June/2022
## Additional Tests
basic_fit = function() {
  load(data_fname)
  fit0 = lm(lvalue ~ cu + ldist + lrgdpcc + lrgdp + regional + comlang + border, data=df)
  return(fit0)
}

#' Replicate the basic results.
replicate_rose_cameron  = function() {
  # Basic results.
  df = load_data()
  ## Rose, etc.
  fit = lm(lvalue ~ cu + ldist + lrgdpcc + lrgdp + regional + comlang + border, data=df)
  summary(fit)  # 1.36 cu
  ## Cameron # Table 3A
  fit = lm(lvalue ~ (cu + ldist + lrgdpcc + lrgdp), data=df)
  summary(fit) # b=1.86
  
  require(sandwich)
  sqrt(vcovCL(fit, cluster=df$ccode1))[2,2]  # 0.4269 -- clustered by CC1
  sqrt(vcovCL(fit, cluster=df$ccode2))[2,2]  # 0.77 -- clustered by CC2
  
  ## without dupes
  load(data_fname)
  fit = lm(lvalue ~ (cu + ldist + lrgdpcc + lrgdp), data=df)
  summary(fit) # b=1.86
  sqrt(vcovCL(fit, cluster=df$ccode1))[2,2]  # 0.2823013 -- clustered by CC1
  sqrt(vcovCL(fit, cluster=df$ccode2))[2,2]  # 0.6235455 -- clustered by CC2
  
  par(mfrow=c(2, 2))
  plot(fit, pch=20, cex=0.3, col=df$comlang+1)
  
  # with dupes removed
  df = load_data(remove_dup = T)
  # Original
  fit0 = lm(lvalue ~ cu + ldist + lrgdpcc + lrgdp + regional + comlang + border, data=df)
  summary(fit0) # 1.06
  AIC(fit0)
  
  fit1 = lm(lvalue ~ (cu + ldist + lrgdpcc + lrgdp + regional + comlang + border)^2, data=df)
  summary(fit1) # 1.06
  AIC(fit1)
}

Q_2way_missing = function() {
  
  df = load_data(remove_dup = T)
  fit = lm(lvalue ~ cu + ldist + lrgdpcc + lrgdp + regional + comlang + border, data=df)
  summary(fit)
  
  ## 1. Is there a two-way structure? 
  
  ## Construct e_ij matrix.
  countries = get_countries(df)
  N = length(countries) # 126 countries.
  E = matrix(NA, nrow=N, ncol=N)
  for(i in 1:nrow(df)) {
    row = df[i, 1:2]
    I = match(row, countries)
    E[I[1], I[2]] = fit$residuals[i]
    E[I[2], I[1]] = fit$residuals[i]
  }
  mean(E, na.rm=T)  # should be 0.
  
  #
  fit1 = lm(fit$residuals ~ as.factor(df$ccode1))
  summary(fit1)
  # Rejects F-test. So, there is some cluster structure.
  
  fit2 = lm(fit1$residuals ~ as.factor(df$ccode2)) #
  summary(fit2)  # F-test rejects. So there is residual cluster structure.
  
  # Yes, there is some cluster structure.
  # 2. missingess
  A = !is.na(E)  # observation pattern.
  library(igraph)
  g = graph_from_adjacency_matrix(A, mode="undirected")
  plot(g, vertex.size=13, vertex.label=countries, 
       edge.arrow.size=0.1, col=as.factor(df$cont1))
  
  imax = which.max(rowSums(A))
  countries[imax]
  
  d = degree(g)
  hist(d, breaks=50)
  countries[tail(order(d), 10)] # most connected countries
  countries[head(order(d), 20)] # least connected countries
  
  # currency union
  in_cu = unique(c(subset(df, cu==1)$ccode1, subset(df, cu==1)$ccode2))
  I_ci = match(in_cu, countries)
  hist(d[I_ci], add=T, col=rgb(1, 0, 0, alpha=0.2), breaks=50)
  
  # Conclusion => Countries in CU are at the left-tail of missingness.
  # i.e., we observe more pairs between countries not in a currency union.
  
  deg = sapply(1:nrow(df), function(i) {
    I = match(c(df$ccode1[i], df$ccode2[i]), countries)
    mean(d[I])
  })
  df1 = cbind(df, deg=deg)
  boxplot(df1$deg ~ df1$cu)
  # deg = average degree of TRADE pair.
  fit1 = lm(lvalue ~ cu + ldist + lrgdpcc + lrgdp + regional + comlang + border + as.factor(deg > 60), data=df1)
  summary(fit1) # CU is not significant anymore!
  
  fit1 = lm(lvalue ~ cu + ldist + lrgdpcc + lrgdp + as.factor(deg> 60), data=df1)
  summary(fit1) # CU is significant again. Common language is significant.
  
  table(df$comlang, df$cu) # common language is important for CU in the dataset.
}

#' Creates M = NxN matrix where countries (i,j) have same colname value.
#' M_ij =1 if i and j have same value (e.g., language)
create_dyadic_filter = function(df, colname) {
  countries = get_countries(df)
  N = length(countries) # 126 countries.
  stopifnot(all(colname %in% colnames(df)))
  stopifnot(all(df[[colname]] %in% c(0, 1))) # should be logical
  
  m = matrix(0, nrow=N, ncol=N)
  for(i in 1:N) {
    j = intersect(which(df$ccode1_i==i), which(df[[colname]]==1))
    match_i = df$ccode2_i[unique(j)]
    m[i, match_i] = 1
    m[match_i, i] = 1
  }
  return(m) # NxN matrix where countries (i,j) have same colname value.
}

CHECK_create_dyadic_filter = function() {
  load(data_fname)
  m = create_dyadic_filter(df, colname = "comcont")
  print(paste(m[3,13], "", m[13,3]))
  countries = get_countries(df)
  countries[3]
  countries[13]
}

#' Transforms a column from <DF> into a NxN matrix (symmetric).
create_dyadic = function(df, colname) {
  N = length(get_countries(df))
  m = matrix(NA, nrow=N, ncol=N)
  stopifnot(colname %in% colnames(df))
  for(i in 1:nrow(df)) {
    I = c(df$ccode1_i[i], df$ccode2_i[i])
    m[I[1], I[2]] = df[[colname]][i]  # Y_ij 
    m[I[2], I[1]] = df[[colname]][i]
  }
  return(m)
}

#' Given dataset, produces decomposition of units = { set_1, set_2, ...} such that for every set_j we observe values for all pairs.
#' filter_by = (optional) vector of (binary) column names to filter on.
#'             e.g., if comm_continent is in the filter, then set_j contains only countries within the same continent.
#'             
create_decomposition = function(df, filter_by=c()) {
  # Represent missing pair data as a graph.
  library(igraph)
  library(nnet)
  Y = create_dyadic(df, "lvalue")
  A = !is.na(Y)  # observation pattern. NxN (country x country) 1=obs pair, 0=not obs
  A0 = A  # 9230 edges
  print(paste("> Total observed country pairs = ", sum(A0)))
  
  # Apply the filter.
  if(length(filter_by) > 0) {
    stopifnot(all(is.character(filter_by)))
    for(f in filter_by) {
      print(paste(">>> Applying Filter by ", f))
      M0 = create_dyadic_filter(df, f)
      A0 = A0 * M0
      print(paste(">>> Total observed country pairs = ", sum(A0)))
    }
  }
  
  # Then find maximum cliques.
  decomposition = list()
  STOP = FALSE
  while(!STOP) {
    g0 = graph_from_adjacency_matrix(A0, mode="undirected")
    out = max_cliques(g0)
    len = unlist(lapply(out, length))
    
    # print(">>> NEW RUN: Distribution of maximal clique length:")
    # print(table(len))
    if(all(len <= 1)) {
      STOP = TRUE
    } else {
      I = as.numeric(out[[which.is.max(len)]]) # break ties at random.
      L = length(decomposition)
      stopifnot(length(intersect(I, unlist(decomposition)))==0) # make sure we are adding disjoint sets
      decomposition[[L + 1]] = I
      A0[I, ] = 0 # remove matched units.
      A0[ ,I] = 0
    }
  }
  # decomposition = LIST()  -- decomposes units into cliques. Assumes exchangeability within cliaues.
  # This can change. We can try different ideas here.
  num = length(unlist(decomposition))
  print(paste(">>> Total cliques=", length(decomposition), "total countries=", num, "(", round(100*num/nrow(Y), 2), "%)"))
  print("> Statistics on clique size:")
  sizes = unlist(lapply(decomposition, length))
  print(summary(sizes))
  
  return(decomposition)
}

#' Checks for create_decom
CHECK_decomposition = function() {
  load(data_fname)
  dec = create_decomposition(df)
  for(j in 1:length(dec)){
    I = dec[[j]]
    E = subset(df, ccode1_i %in% I  & ccode2_i %in% I) # 
    # all pairs in I should be observed.
    stopifnot(nrow(E) == length(I)*(length(I)-1)/2) ## must be a clique.
    print("  [ PASS ]  ")
  }
  
  load(data_fname)
  dec = create_decomposition(df, filter_by = c("comcont", "comlang"))
  for(j in 1:length(dec)){
    I = dec[[j]]
    E = subset(df, ccode1_i %in% I  & ccode2_i %in% I) # 
    # all pairs in I should be observed.
    stopifnot(nrow(E) == length(I)*(length(I)-1)/2) ## must be a clique.
    stopifnot(all(E$comcont==1))
    stopifnot(all(E$comlang==1))
    print("  [ PASS ]  ")
  }
  
  
}

#' Centers a matrix that may have NA values (missing data.)
#' 
double_center = function(M, nreps=200) {
  
  M = M - mean(M, na.rm=T)
  return(M)
  # Old code. Not necessary to center in rows-columns.
  # once = function(A) {
  #   Ind = is.na(A)
  #   A[Ind] = 0
  #   N = nrow(A)
  #   stopifnot(nrow(A) == ncol(A))
  #   #U = matrix(1, nrow=N, ncol=N)
  #   # A - U %*% A/N - A %*% U/N + mean(A)*U
  #   r = rowMeans(A)
  #   c = colMeans(A)
  #   ones = rep(1, N)
  #   U = ones %*% t(ones)
  #   
  #   ret = A - r %*% t(ones) - ones %*% t(c) + mean(A) * U
  #   ret[Ind] = NA
  #   return(ret)
  # }
  # 
  # B = M
  # for(j in 1:nreps) {
  #   B = once(B)
  # }
  # 
  # lmin = max(abs(c(rowMeans(B, na.rm=T), colMeans(B, na.rm=T))))
  # print(paste("Centered for ", nreps, "iterations. ", "Max row/col mean=", lmin))
  # return(B)
}

# transform between matrix and vector format.
# to_v(to_m(A)) = A
to_v = function(m) m[lower.tri(m)]
to_m = function(v, N) {
  cf = c(-2*length(v), -1, 1); 
  K = polyroot(cf)
  K = Re(K)[which(Re(K) > 0)]
  if(K != N) {warning("Not using N=126 as country index")}
  m = matrix(NA, nrow=K, ncol=K); 
  m[lower.tri(m)] = v; 
  m1 = t(m)
  m[upper.tri(m)] = m1[upper.tri(m1)]
  return(m)
}

#' Transforms the data into the appropriate format.
#' 
preprocess_RRI = function(beta_0, center,filter_cols=c()) {
  
  # Load data (df)
  load(file=data_fname)

  # Get countries
  countries = get_countries(df)
  stopifnot(max(c(df$ccode1_i, df$ccode2_i))==length(countries))
  N = length(countries) # 126 countries.
  
  # Transform data to dyadic.
  #' 
  # Y = A X_cu + B X_controls    # nx1 , nx
  Y = create_dyadic(df, "lvalue")
  stopifnot(all(Y==to_m(to_v(Y), N), na.rm=T))
  
  CU = create_dyadic(df, "cu")
  # controls
  LDIST = create_dyadic(df, "ldist")
  LRGDP = create_dyadic(df, "lrgdp")
  LRGDPCC = create_dyadic(df, "lrgdpcc")
  COMLANG = create_dyadic(df, "comlang")
  REGIONAL = create_dyadic(df, "regional")
  
  # center?
  if(center) {
    print(" Centering...")
    Y = double_center(Y)
    CU = double_center(CU)
    LDIST  = double_center(LDIST )
    LRGDP = double_center(LRGDP)
    LRGDPCC = double_center(LRGDPCC)
    COMLANG  = double_center(COMLANG )
    REGIONAL  = double_center(REGIONAL )
  }
  
  # Make into a long vector.
  y = to_v(Y)
  # currency union index
  X_cu = to_v(CU)
  # controls
  X_controls = cbind(ldist=to_v(LDIST), 
                     lrgdp=to_v(LRGDP),
                     lrgdpcc=to_v(LRGDPCC),
                     comlang=to_v(COMLANG),
                     regional=to_v(REGIONAL))
  
  ## Restricted model (Assume Beta_CU=0)
  fit0 = lm(y ~ X_controls)
  
  # Create E = NxN matrix of residuals. Need to be careful about NAs.
  ObsIndex = to_v(!is.na(Y)) # length = N(N-1)/2  1=number, 0=NA
  E = rep(NA, length(ObsIndex))
  stopifnot(length(fit0$residuals)==sum(ObsIndex))
  E[ObsIndex] = fit0$residuals
  E = to_m(E, N)            # Contains residuals in NxN form. E_ij not NA if (i, j) observed.
  ## E[3,121]
  ## CHECKING E: All observed data should correspond to numerical values in E.
  for(irow in 1:nrow(df)) {
    stopifnot(!is.na(E[df$ccode1_i[irow], df$ccode2_i[irow]])) 
  }
  stopifnot(all(to_v(!is.na(Y)) == to_v(!is.na(E))))  # 
  
  # Generate decomposition of units = {set_1, set_2, ...} such that for every set_j we observe values for all pairs.
  print(paste("> Generating decomposition..."))
  decomposition = create_decomposition(df, filter_by=filter_cols)
  
  list(E=E, y=y, X_cu=X_cu, X_controls=X_controls, dec=decomposition)
}

#' Tests H0: beta_CU=..  via residual randomization.
RRI_trade = function(beta_0, obj_preprocess, num_r=2000) {
  
  E = obj_preprocess$E # residual matrix (incl. NA)
  y = obj_preprocess$y # outcome vector
  X_cu = obj_preprocess$X_cu # CU indicator 
  X_controls = obj_preprocess$X_controls # controls 
  decomposition = obj_preprocess$dec # decomposition of units into "cliques" (no NA for every pair in a clique.)
  
  # OLS 
  fit = lm(y ~ X_cu + X_controls)
  summary(fit)  # exacty the same as Cameron
  Tn = coef(fit)[2] - beta_0 # test statistic
  
  na.index = is.na(E)
  
  tvals = replicate(num_r, {
    E_r = E
    for(j in 1:length(decomposition)) {
      I_j  = decomposition[[j]]
      stopifnot(sum(is.na(E_r[I_j, I_j])) > 0) # only permute non-NA values.
      ## Shuffle
      I_j_new = sample(I_j) 
      E_r[I_j, I_j] = E_r[I_j_new , I_j_new] # shuffle rows and columns.
    }
    # Check if we "respected" the NA structure.
    na.index_r = is.na(E_r)
    stopifnot("E and E_r should have the same NAs"=(all(na.index_r==na.index)))
    
    fit_r = lm(to_v(E_r) ~ X_cu  + X_controls)
    
    coef(fit_r)[2]
  })
  
  # hist(tvals, breaks=50, main=paste("RRI double perm: b^=", round(coef(fit)[2],3), "- beta_0=", beta_0))
  # abline(v=tobs, col="red")
  randomization_pval(Tn, tvals, two.sided = T)
}


RRI_trade_no_structure = function(beta_0, obj_preprocess, type, num_r=5000) {
  # type = {perm, sign}
  E = obj_preprocess$E
  y = obj_preprocess$y
  X_cu = obj_preprocess$X_cu
  X_controls = obj_preprocess$X_controls
  decomposition = obj_preprocess$dec
  
  fit = lm(y ~ X_cu + X_controls)
  summary(fit)  # exacty the same as Cameron
  tobs = coef(fit)[2] - beta_0
  
  e_r = to_v(E)
  g = create_invariance(length(e_r), type=type) 
    
  tvals = replicate(num_r, {
    # E_r residuals.
    fit_r = lm(g(e_r) ~ X_cu  + X_controls)
    coef(fit_r)[2]
  })
  
  randomization_pval(tobs, tvals, two.sided = T)
}

# Single simulation
single_sim = function(beta1_H0, num_r, center, filters) {
  
  out = preprocess_RRI(center = center, filter_cols=filters)
  rri_test = function(b) { 
    print(paste(">>> Testing H0: beta=", b))
    RRI_trade(beta_0=b, out, num_r=num_r) 
  }
  
  rri_test(beta1_H0) # p.values
}

# Main simulation
main_sim = function(beta1_H0, num_r) {
  
  cols = c("center", "same_cont", "same_lang", "same_border", "dec_size", "dec_mean", "dec_sd", "beta1_H0", "p.value")
  Results = matrix(0, nrow=0, ncol=length(cols))
  colnames(Results) = cols
  
  center_vals = c(T) # always center
  filter_vals = list(c(), c("comcont"), c("comlang"), c("border"), 
                     c("comcont", "comlang"), c("comcont", "border"), c("comlang", "border"))
                     
  for(center in center_vals) {
    for(ii in 1:length(filter_vals)) {
      
      f = filter_vals[[ii]]
      print(paste("Center? = ", center))
      print("Filters = ")
      print(f)
      ## 
      out = preprocess_RRI(center = center, filter_cols=f)
      pval = single_sim(beta1_H0, num_r, center, f)
      
      dec_vals = unlist(out$dec)
      
      Results = rbind(Results, as.numeric(c(center, "comcont" %in% f, "comlang" %in% f, "border" %in% f,
                                            length(dec_vals), round(mean(dec_vals), 2),  round(sd(dec_vals), 2),
                                            beta1_H0, round(pval, 3))))
      
      print("-==-=-=-=-=-===-   RESULTS   -=-=-=- -=-= -=- =-")
      print(Results)
    }
    
  }
  
  # Store results?
  return(as.data.frame(Results))
}

# Main simulation with simple structure
main_sim_no_dyadic = function(num_r) {
  set.seed(539)
  types = c("perm", "sign", "double")
  for(t in types) {
    print("-=-=-=-==-=-=-=-==-=-=")
    print(paste("Testin' invariance type =", t))
    print("-=-=-=-==-=-=-=-==-=-=")
    
    out = preprocess_RRI(center = T, filter_cols=c())
    H0 = function(b) {
      RRI_trade_no_structure(b, out, type = t, num_r=num_r)
    }
    
    invert_rrtest(H0, CI_initial = c(-3, 3))
  }
  
}

# some quick examples
quick_examples = function() {
  set.seed(555)
  out = preprocess_RRI(center = T, filter_cols=c("comcont"))
  # simple exchangeability rejects Beta_CU=0
  RRI_trade_no_structure(beta_0=0, out, type = "perm", num_r=3000)
  
  # dyadic exchangeability does not reject
  RRI_trade(beta_0=0, out, num_r = 3000)
}

analyze_results_From_mercury = function() {
  load("CI_5.3.rda")
  L = 1 # determines how we accept a null value.
       # smaller means we are more conservative
  # Only results for no-filter model are affected.
  
  for(j in 1:length(ci)) {
    m1 = ci[[j]]
    model = loess(pval ~ theta0, m1, span=0.05)
    
    plot(m1$theta0, m1$pval, pch=20, main=paste("Scenario", j))
    lines(m1$theta0, predict(model, m1$theta0), lwd=2, col="blue")  
    
    Acc = (m1$pval >= 0.05)
    Acc_near = rep(1, length(Acc))
    if(L > 0) {
      Acc_near = sapply(seq(1+L, length(m1$pval)-L), function(i) {
        Nei = c(seq(i-L, i-1), seq(i+1, i+L))
        prod(Acc[Nei])
      })
      Acc_near = c(rep(0, L), Acc_near, rep(0, L))
    }
    Ij = which(Acc*Acc_near==1)
    left = round(min(m1$theta0[Ij]), 4)
    right = round(max(m1$theta0[Ij]), 4)
    
    print(paste("CI for scenario ", j, " = [", left, ",", right, "]"))
  }
  
  
}

#' Some additional tests on the sensitivity of results.
#' The code gives some additional evidence why the coefficient of currency union is probably not 0.
sensitivity_analysis = function(){
  
  load(data_fname)
  
  M = !is.na(create_dyadic(df, "lvalue")) # missingness. 9230 observed.
  print(paste("observed", sum(M), "out of ", prod(dim(M)))) # 58%
  countries = get_countries(df)
  rowSums(M) # degree of every country
  df$deg = 0.5*(rowSums(M)[df$ccode1_i] + rowSums(M)[df$ccode2_i]) # average degree of dyad
  boxplot(df$deg ~ df$cu) # Low connectivity almost perfect predictor of CU=1
  
  df = subset(df, ccode1 != "CMR" & ccode2 != "CMR") # 50% of CU=1 is CMR.
  # CMR    GNQ
  fit0 = lm(lvalue ~ cu + ldist + lrgdpcc + lrgdp + regional + comlang + border, data=df)
  summary(fit0)
  # effect goes away with weighted LS
  fit1 = lm(lvalue ~ cu + ldist + lrgdpcc + lrgdp + regional + comlang + border, data=df, weights=1/(0.2+df$cu))
  summary(fit1)
  AIC(fit1)
  tobs = coef(fit0)[2]
  
  # FRT based on prop scores
  fit_ps = glm(cu ~ ldist + lrgdpcc + lrgdp + regional + comlang + border, data=df, family="binomial")
  phat = fit_ps$fitted.values
  tvals = replicate(2000, {
    cu1 = rbinom(nrow(df), size=1, prob=phat)
    fit1 = lm(lvalue ~ cu1 + ldist + lrgdpcc + lrgdp + regional + comlang + border, data=df)
    coef(fit1)[2]
  })
  par(mfrow=c(1, 1))
  hist(tvals)
  abline(v=tobs, col="red")
  randomization_pval(tobs, tvals, T) # don't reject that beta_CU=0
  
  # DiCiccio and Romano Test based on permutations.
  tstat = function(lm_obj) {
    val = coef(lm_obj)[2]
    se = summary(lm_obj)$coefficients[2,2]
    val/se
  }
  tobs = tstat(fit0)
  
  tvals = replicate(2000, {
    cu1 = sample(df$cu)
    fit1 = lm(lvalue ~ cu1 + ldist + lrgdpcc + lrgdp + regional + comlang + border, data=df)
    tstat(fit1)
  })
  hist(tvals, breaks=100)
  abline(v=tobs, col="red")
  randomization_pval(tobs, tvals, T) # reject with pval < 0.04 (weak evidence.)
  
  # DiCiccio and Romano Test based on residuals
  tobs = tstat(fit0)
  e = fit0$residuals
  I1 = which(df$cu==1)
  I0 = which(df$cu==0)
  summary(e[I1])
  summary(e[I0])
  hist(e[I1])
  hist(e[I0])
  qqplot(e[I0], e[I1])
  abline(0, 1, col="red")
  # extrema in the residuals do not line up. ( between e|cu=1  and e|cu=0  )
  
  # fit under the null
  fit0 = lm(lvalue ~ ldist + lrgdpcc + lrgdp + regional + comlang + border, data=df)
  yhat = fit0$fitted.values
  #
  fit1 = lm(lvalue ~ cu + ldist + lrgdpcc + lrgdp + regional + comlang + border, data=df)
  tobs = tstat(fit1)
  tvals = replicate(2000, {
    y1 = yhat + sign(rnorm(length(y0)))*e
    fit1 = lm(y1 ~ cu + ldist + lrgdpcc + lrgdp + regional + comlang + border, data=df)
    tstat(fit1)
  })
  hist(tvals, breaks=100)
  abline(v=tobs, col="red")
  randomization_pval(tobs, tvals, T) # does not reject.
  
  
  ## Leverage of design
  fit0 = lm(lvalue ~ cu + ldist + lrgdpcc + lrgdp + regional + comlang + border, data=df)
  
  max(hatvalues(fit0)) / mean(hatvalues(fit0))
  
  # continent is predictive of residuals.
  par(mfrow=c(2, 2))
  plot(fit0)
  e = fit0$residuals
  summary(lm(e ~ factor(df$cont1)))
  
  boxplot(e ~ df$cont1)
  
  #' In summary
  #'  - weighted LS does not reject beta_CU=0
  #'  - CU is nearly perfectly correlated with trade degree.
  #'  - High leverage (only 21 datapoints with CU=1 out of thousands)
  #'  - FRT-based randomization does not reject.
  #'  - residual permutation test a-la Romano does not reject.
  #'  - permutation test a-la Romano rejects with p-val=0.04
}
