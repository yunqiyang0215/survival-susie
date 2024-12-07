# This script is a first draft of the CoxPH-SuSiE demo.
# remotes::install_github("karltayeb/logisticsusie")
# remotes::install_github("stephenslab/susieR")
library(survival)
library(susieR)
library(logisticsusie)

# Function to calculate approximate log-Bayes factor.
# z: zscore of the regression coefficient
# s: s.d. of the estimated coefficient.
# v0: prior variance
compute_lbf <- function (z, s, v0)
  log(s^2/(s^2 + v0))/2 + z^2/2*v0/(s^2 + v0)

compute_approx_post_var <- function (z, s, v0)
  1/(1/s^2 + 1/v0)

# @param v1: posterior variance
# @param s: s.d. of the estimated coefficient
# @param bhat: estimated beta effect
compute_approx_post_mean <- function (v1, s, bhat)
  v1*bhat/s^2

# TO DO: Improve the code for this function a bit.
surv_uni_fun <- function (x, y, o, v0, estimate_intercept = 0, ...) {
  fit  <- coxph(y ~ x + offset(o))
  bhat <- summary(fit)$coefficients[1,1]
  sd   <- summary(fit)$coefficients[1,3]
  z    <- bhat/sd
  lbf  <- compute_lbf(z,sd,v0)
  lbf_corr <- lbf - z^2/2 + summary(fit)$logtest[1]/2
  var  <- compute_approx_post_var(z,sd,v0)
  mu   <- compute_approx_post_mean(var,sd,bhat)
  return(list(mu = mu,var = var,lbf = lbf_corr,intercept = 0))
}

L <- 3

load("data/survival_demo.RData")

# Create the survival object.
pheno <- Surv(time,status)

fit_coxph <- ser_from_univariate(surv_uni_fun)

# This step may take a minutes or two to run
fit <- ibss_from_ser(geno,pheno,L = L,ser_function = fit_coxph,
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
