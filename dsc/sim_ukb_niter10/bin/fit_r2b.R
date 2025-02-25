
library(R2BGLiMS)


p = ncol(dat$X)
dat2 = as.data.frame(dat$X)
predictors = colnames(dat$X)
dat2$surT = dat$y[,1]/max(dat$y[,1])
dat2$status = dat$y[,2]

fit <- R2BGLiMS(
  likelihood="Weibull",
  data=dat2,
  outcome.var="status",
  times.var="surT",
  model.space.priors=list(list("a"=1,"b"=length(predictors),"Variables"=predictors))# Beta-binomial(1,P) model space prior
)

pip <- fit@posterior.summary.table[3:(2+p), 'PostProb']
effect_estimate <- fit@posterior.summary.table[3:(2+p), 'Mean']

