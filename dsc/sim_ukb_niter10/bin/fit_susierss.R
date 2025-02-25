# load packages
# Modified Karl's code for intercept part
library(survival)
library(susieR)
devtools::load_all("/project2/mstephens/yunqiyang/surv-susie/logisticsusie")


get_coxph_sumstats = function(X, y){
  p = ncol(X)
  bhat = rep(NA, p)
  sebhat = rep(NA, p)
  for (j in 1:p){
    fit = coxph(y ~ X[ ,j])
    bhat[j] = summary(fit)$coefficients[1, 1]
    sebhat[j] = summary(fit)$coefficients[1, 3]
   #zscore[j] = summary(fit)$coefficients[1, 4] 
 }
 return(list(bhat = bhat, sebhat = sebhat))
} 



y <- Surv(dat$y[,1], dat$y[,2])
X = as.matrix(dat$X)

n = nrow(X)
R = cor(X)
sumstats = get_coxph_sumstats(X, y)
niter = 100


fit <- susie_rss(bhat = sumstats$bhat, shat = sumstats$sebhat, n = n, R = R, var_y = 1, L = 5,
                         estimate_residual_variance = FALSE, niter = niter)

pip <- logisticsusie:::get_pip(fit$alpha)
effect_estimate <- colSums(fit$alpha * fit$mu)
class(fit) = "susie"
cs <- susie_get_cs(fit, X)
iter <- fit$niter

