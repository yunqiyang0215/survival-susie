surv_uni_fun <- function (x, y, e, v0, estimate_intercept = 0, ...) {
  fit  <- coxph(y ~ x + offset(e))
  out  <- summary(fit)$coefficients
  bhat <- out[1,"coef"]
  s    <- out[1,"se(coef)"]
  z    <- bhat/s
  lbf  <- log(s^2/(v0 + s^2))/2 + z^2/2*v0/(v0 + s^2)
  lbf  <- lbf - z^2/2 + summary(fit)$logtest[1]/2
  v1   <- 1/(1/v0 + 1/s^2)
  mu1  <- v1*bhat/s^2
  return(list(mu = mu1,var = var1,lbf = lbf,intercept = 0))
}

fit_coxph <- ser_from_univariate(surv_uni_fun)

# This step may take a minutes or two to run
fit <- ibss_from_ser(geno,pheno,L = 3,ser_function = fit_coxph,
                     tol = 0.0001,maxit = 100)

class(fit) <- c("susie","list")
out <- susie_get_cs(fit,geno,min_abs_cor = 0.1)

# CS 1.
i <- 555
print(b[i])
print(fit$mu[1,i])

# CS 2.
i <- 994
print(b[i])
print(fit$mu[2,i])

# TO DO:
#
# (1) Create a plot showing the coxph-based association p-values.
# (1) Create a nice PIP plot showing the CSs and the causal SNPs.
# (2) Create a scatterplot of true vs. estimated coefs.
# 
