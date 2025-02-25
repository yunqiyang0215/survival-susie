library(survival.svb)

X = as.matrix(dat$X)
y = dat$y[,1]
delta = dat$y[,2]
fit <- survival.svb::svb.fit(y, delta, X, maxiter = niter)

pip <- fit$inclusion_prob
effect_estimate <- fit$beta_hat
