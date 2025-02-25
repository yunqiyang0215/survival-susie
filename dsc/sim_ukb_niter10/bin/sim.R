library(Matrix)
library(survival)


sim_surv <- function(b, X, censor_lvl, dist = "exp"){
  n = nrow(X)
  p = ncol(X)
  
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
  
  dat = list()
  dat[["X"]] = X
  dat[["y"]] = y
  return(dat)
}


geno <- readRDS("/project2/mstephens/yunqiyang/surv-susie/realdata/bloodcells_chr3.18510058.19066126.rds")

# Select a random sample of 50k.
n = 5e4
indx = sample(1:nrow(geno), n, replace = FALSE)
geno_sub = geno[indx, ]


# Remove SNPs with MAF < 0.01
maf = apply(geno_sub, 2, function(x) sum(x)/2/length(x))
geno_sub2 = geno_sub[, maf > 0.01]
snp_total = ncol(geno_sub2)

# Start from a random point on the genome
indx_start = sample(1: (snp_total - p), 1)
X = geno_sub2[, indx_start:(indx_start + p -1)]

# create effect size vector b.
b = rep(0, p)
is_effect = rep(0, p)
indx_snp = 0

if (num_effect != 0) {
  indx_snp <- sample(2:ncol(X), num_effect, replace = TRUE)
  if (cor_type == "independent"){
    # permute X to make them nearly independent
    for (i in 1:p){
      X[, i] = X[, i][sample(1:n, replace = FALSE)]
    }
  }
  effect_size = rnorm(num_effect, mean = 0, sd = sqrt(prior_variance))
  b[indx_snp] = effect_size
  is_effect[indx_snp] = 1
}


dat = sim_surv(b = c(1, b), X, censor_lvl)
X = dat$X
y = dat$y










