
#### Helper functions ######
# Function to calculate approximate BF based on Wakefield approximation
# @param z: zscore of the regression coefficient
# @param s: standard deviation of the estimated coefficient
compute_abf <- function(z, s, prior_variance){
  abf <- sqrt(s^2/(s^2+prior_variance))*exp(z^2/2*(prior_variance/(s^2+prior_variance)))
  return(abf)
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
  bf <- compute_abf(zscore, sd, prior_variance)
  var <- compute_approx_post_var(zscore, sd, prior_variance)
  mu <- compute_approx_post_mean(var, sd, bhat)
  lbf <- log(bf)
  return(list(mu = mu, var=var, lbf=lbf, intercept=0))
}

fit_coxph <- ser_from_univariate(surv_uni_fun)
