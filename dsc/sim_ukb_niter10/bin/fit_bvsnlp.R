library(BVSNLP)
library(survival)

# Create survival response
y <- Surv(dat$y[,1], dat$y[,2])
X = as.data.frame(dat$X)
p = ncol(X)

fit = bvs(X, resp = y, prep = FALSE, family = "survival", mod_prior = "beta")

effect_estimate = rep(0, p)
effect_estimate[fit$HPM] = fit$beta_hat
pip = fit$inc_probs

