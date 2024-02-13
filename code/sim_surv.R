# Function to simulate survival time under exponential model. That is,
# assuming survival time is exponentially distributed.
# lambda(t) = lambda*exp(b0 + Xb). S(t) = exp(-lambda*t),
# F(t) = 1- S(t) \sim Unif(0,1). Therefore, t = log(1-runif(0,1))/-exp(b0+Xb).
# For censored objects, we simulate the censoring time by rescale the actual survival time.
# @param b: vector of length (p+1) for true effect size, include intercept.
# @param X: variable matrix of size n by p.
# @param censor_lvl: a constant from [0,1], indicating the censoring level in the data.
# @return  dat: a dataframe that contains `y`, `x` and `status`.
# `status`: censoring status: 0 = censored, 1 = event observed. See Surv() in library(survival)
sim_surv_exp <- function(b, X, censor_lvl){
  n = nrow(X)
  p = ncol(X)
  dat = list()

  status <- ifelse(runif(n) > censor_lvl, 1, 0)
  lambda <- exp(cbind(rep(1,n), X) %*% b)
  surT <- log(1 - runif(n)) /(-lambda)
  # rescale censored subject to get observed time
  surT[status == 0] = surT[status == 0] * runif(sum(status == 0))

  y = cbind(surT, status)
  colnames(y) = c("time", "status")
  colnames(X) <- unlist(lapply(1:p, function(i) paste0("x", i)))
  dat[["X"]] = X
  dat[["y"]] = y
  return(dat)
}
