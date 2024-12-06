## censor_lvl: 0, 0.2, 0.4, 0.6, 0.8

# Get the genotype data.
X <- readRDS("../data/Thyroid_ENSG00000132855.rds")$X
rownames(X) <- NULL
colnames(X) <- NULL

# Remove SNPs with minor allele frequencies less than 1%.
maf <- colMeans(X)/2
maf <- pmin(maf,1 - maf)
X   <- X[,maf > 0.01]

# Keep 1,000 of the SNPs.
X <- X[,3501:4500]

# Simulate the coefficients, b.

sim_surv <- function(b, X, censor_lvl, dist = "exp") {
  n <- nrow(X)
  p <- ncol(X)
  
  # survival distribution with rate lambda.surv
  lamb.surv <- exp(cbind(rep(1,n), X) %*% b)
  lamb.mean <- mean(lamb.surv)
  
  # censoring distribution with rate lambda.censor
  lamb.censor <- censor_lvl*lamb.mean/(1-censor_lvl)
  
  # generate survival time and censor time
  survT <- log(1 - runif(n)) /(-lamb.surv)
  censorT <- log(1 - runif(n)) /(-lamb.censor)
  
  status <- ifelse(survT < censorT, 1, 0)
  time = pmin(survT, censorT)
  y = cbind(time, status)
  colnames(y) = c("time", "status")
  colnames(X) <- unlist(lapply(1:p, function(i) paste0("x", i)))
}


# create effect size vector b.
b = rep(0, p)
is_effect = rep(0, p)
indx_snp = 0

if (num_effect != 0) {
  if (cor_type == "real"){
    indx_snp <- sample(2:n, num_effect, replace = FALSE)
  }
  effect_size = rnorm(n = num_effect)
  b[indx_snp] = effect_size
  is_effect[indx_snp] = 1
}

dat = sim_surv(b = c(1, b), X, censor_lvl)

