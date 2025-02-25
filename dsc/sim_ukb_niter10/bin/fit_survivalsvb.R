library(survival.svb)

X = as.matrix(dat$X)
y = dat$y[,1]
delta = dat$y[,2]
maxiter = 100
fit <- survival.svb::svb.fit(y, delta, X, maxiter = maxiter)

pip <- fit$inclusion_prob
effect_estimate <- fit$beta_hat
