
### Helper functions to compute BF
compute_multiple_lbf <- function(x, y, o, prior_variance, ...){
  fit <- coxph(y~ x + offset(o))
  bhat <- summary(fit)$coefficients[1, 1] # bhat = -alphahat
  sd <- summary(fit)$coefficients[1, 3]
  zscore <- bhat/sd
  t1 = proc.time()
  lbf <- compute_lbf(zscore, sd, prior_variance)
  lbf.corr <- lbf - bhat^2/sd^2/2+ as.numeric(summary(fit)$logtest[1]/2)
  t2 = proc.time()
  return(list(zscore = zscore, bhat = bhat, sd = sd,
              lbf=lbf, lbf.corr = lbf.corr, time.lp = (t2 - t1)[3]))
}

### Numerical integration

#' Gaussian quadrature
#'
#' approximate log(\int f(x)dx) by quadrature on \int f(x)/p(x) p(x)
#' where p(x) is a normal density, we approximate this integral
#' with a gaussian quadrature rule via `statmod::gauss.quad.prob`
#' done on a log scale for numerical stability (f is positive)
#' the goal is to choose p(x) such that f(x)/p(x) is not too variable near the mode
#' @param mu mean of gaussian quadrature
#' @param sd sd of gaussian quadrature
gauss_quad2 <- function(log_f, mu, sd, n=32){
  q <- statmod::gauss.quad.prob(n, dist='normal', mu=mu, sigma=sd)
  log_integrand <- function(x){log_f(x) - dnorm(x, mu, sd, log=T)}
  return(matrixStats::logSumExp(log_integrand(q$nodes) + log(q$weights)))
}

# @param b: the effect size
# @param surT: a Surv() object, containing time and status
# @param x: covariate
get_partial_loglik <- function(b, survT, x){
  ll.partial = logLik(coxph(survT ~ x, init = c(b), control=list('iter.max'= 0, timefix=FALSE)))
  return(as.numeric(ll.partial))
}

# @param b: the effect size
# @param surT: a Surv() object, containing time and status
# @param x: covariate
# @param sig0: prior standard deviation
get_log_f <- function(b, survT, x, sig0){
  partial_ll = get_partial_loglik(b, survT, x)
  return(partial_ll + dnorm(b, 0, sd = sig0, log = TRUE))
}


########

library(survival)
library(statmod)
source("/project2/mstephens/yunqiyang/surv-susie/survival-susie/code/surv_susie_helper.R")

survT <- Surv(dat$y[,1], dat$y[,2])
n = length(dat$y[,1])
prior_sd = 1


if (flip == "flip"){
  if (binary == "bin"){
    x = 1 - dat$X[,1]
  }else{
    x = 2 - dat$X[,1]
  }
}else{
  x = dat$X[,1]
}

res = surv_uni_fun(x = x, y = survT, o = rep(0, n), prior_variance = prior_sd^2)
mu = res$mu
sd = sqrt(res$var)
log_f <- function(b) get_log_f(b, survT, x = x, sig0 = prior_sd)
t1 = proc.time()
lbf.numerical = gauss_quad2(Vectorize(log_f), mu, sd, n=32) - get_partial_loglik(b = 0, survT, x)
t2 = proc.time()

# compute other lbf
time.gh = (t2 - t1)[3]
res.lbf = compute_multiple_lbf(x = x, y = survT, o = rep(0, n), prior_variance = prior_sd^2)
result = c(res.lbf$zscore, res.lbf$bhat, res.lbf$sd, res.lbf$lbf, res.lbf$lbf.corr, 
           lbf.numerical, res.lbf$time.lp, time.gh)
names(result) = c("zscore", "bhat", "sd", "lbf.wakefeld", "lbf.laplace", "lbf.numerical", "time.lp", "time.gh")

