sim_surv <- function(b, X, censor_lvl, dist = "exp") {
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

####### parameters setting ######
gtex = readRDS("/project2/mstephens/SuSiE/gtex_data/Toys/Thyroid.ENSG00000132855.RDS")

# Remove SNPs with MAF < 0.01
maf = apply(gtex$X, 2, function(x) sum(x)/2/length(x))
X0 = gtex$X[, maf > 0.01]

snp_total = ncol(X0)
n = nrow(X0)

# Start from a random point on the genome
indx_start = sample(1: (snp_total - p), 1)
X = X0[, indx_start:(indx_start + p -1)]

# create effect size vector b.
b = rep(0, p)
is_effect = rep(0, p)
indx_snp = 0

if (num_effect != 0) {
  if (cor_type == "real"){
    indx_snp <- sample(2:n, num_effect, replace = FALSE)
  }
  if (cor_type == "independent"){
    # permute X to make them nearly independent
    for (i in 1:p){
      X[, i] = X[, i][sample(1:n, replace = FALSE)]
    }
    indx_snp <- round(seq(2, p, length.out = num_effect))
  }
  effect_size = rnorm(n = num_effect, mean = 0, sd = sqrt(prior_variance))
  b[indx_snp] = effect_size
  is_effect[indx_snp] = 1
}

dat = sim_surv(b = c(1, b), X, censor_lvl)

