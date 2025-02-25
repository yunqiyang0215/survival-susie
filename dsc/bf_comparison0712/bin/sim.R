sim_surv <- function(b, X, censor_lvl, dist = "exp"){
  n = nrow(X)
  p = ncol(X)
  
  # survival distribution with rate lambda.surv
  lamb.surv <- exp(cbind(rep(1,n), X) %*% b)
  lamb.mean <- mean(lamb.surv)
  
  # censoring distribution with rate lambda.censor
  lamb.censor <- censor_lvl*lamb.mean/(1-censor_lvl)
  
  # generate survival time and censor time
  survT <- log(1 - runif(n)) /(-lamb.surv)
  censorT <- log(1 - runif(n)) /(-lamb.censor)
  
  status <- ifelse(survT < censorT, 1, 0)
  time = pmin(survT, censorT)
  y = cbind(time, status)
  colnames(y) = c("time", "status")
  colnames(X) <- unlist(lapply(1:p, function(i) paste0("x", i)))
  
  dat = list()
  dat[["X"]] = X
  dat[["y"]] = y
  return(dat)
}



# the first element of b is for intercept
b = c(1, b1)
X = cbind(rbinom(n, size = 2, prob = allele_freq))
dat <- sim_surv(b, X, censor_lvl)
