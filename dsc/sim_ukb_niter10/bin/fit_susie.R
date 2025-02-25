# load packages
# Modified Karl's code for intercept part
library(survival)
library(susieR)
devtools::load_all("/project2/mstephens/yunqiyang/surv-susie/logisticsusie")

#### Helper functions ######
# Function to calculate log of approximate BF based on Wakefield approximation
# @param z: zscore of the regression coefficient
# @param s: standard deviation of the estimated coefficient
compute_lbf <- function(z, s, prior_variance){
  abf <- sqrt(s^2/(s^2+prior_variance))
  lbf <- log(sqrt(s^2/(s^2+prior_variance))) + z^2/2*(prior_variance/(s^2+prior_variance))
  return(lbf)
}


compute_approx_post_var <- function(z, s, prior_variance){
  post_var <- 1/(1/s^2 + 1/prior_variance)
  return(post_var)
}

# @param post_var: posterior variance
# @param s: standard deviation of the estimated coefficient
# @param bhat: estimated beta effect
compute_approx_post_mean <- function(post_var, s, bhat){
  mu <- post_var/(s^2)*bhat
  return(mu)
}


surv_uni_fun <- function(x, y, o, prior_variance, estimate_intercept = 0, ...){
  fit <- coxph(y~ x + offset(o))
  bhat <- summary(fit)$coefficients[1, 1] # bhat = -alphahat
  sd <- summary(fit)$coefficients[1, 3]
  zscore <- bhat/sd
  lbf <- compute_lbf(zscore, sd, prior_variance)
  lbf.corr <- lbf - bhat^2/sd^2/2+ summary(fit)$logtest[1]/2
  var <- compute_approx_post_var(zscore, sd, prior_variance)
  mu <- compute_approx_post_mean(var, sd, bhat)
  return(list(mu = mu, var=var, lbf=lbf.corr, intercept=0))
}

fit_coxph <- ser_from_univariate(surv_uni_fun)

####### load data to fit survival susie
#### parameter settings
L = 5
p = ncol(dat$X)
maxiter = 10


## Create  survival object
y <- Surv(dat$y[,1], dat$y[,2])
X = as.matrix(dat$X)


fit <- ibss_from_ser(X, y, L = L, prior_variance = 1., prior_weights = rep(1/p, p), tol = 1e-3, maxit = maxiter,
                     estimate_intercept = TRUE, ser_function = fit_coxph)


pip <- logisticsusie:::get_pip(fit$alpha)
effect_estimate <- colSums(fit$alpha * fit$mu)
class(fit) = "susie"
cs <- susie_get_cs(fit, X)
iter = fit$iter



