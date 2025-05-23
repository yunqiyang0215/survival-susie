# Generate the data set for the CoxPH-SuSiE tutorial.
library(tools)
censor_level <- 0.3

# Set the seed so the results can be reproduced.
set.seed(1)

# Get the genotype data.
geno <- readRDS("../data/Thyroid_ENSG00000132855.rds")

# Remove SNPs with minor allele frequencies less than 1%.
maf  <- colMeans(geno)/2
maf  <- pmin(maf,1 - maf)
geno <- geno[,maf > 0.01]

# Keep 1,000 of the SNPs.
geno <- geno[,3501:4500]

# Simulate the coefficients, b.
n <- nrow(geno)
p <- ncol(geno)
b <- rep(0,p)
i <- c(384,555,994)
b0 <- 0
b[i] <- c(1,-2,2) 
  
# Generate survival times and censor times.
lsurv   <- exp(b0 + drop(geno %*% b))
lcensor <- censor_level * mean(lsurv)/(1 - censor_level)
survT   <- (-log(1 - runif(n))/lsurv)
censorT <- (-log(1 - runif(n))/lcensor)
status  <- as.numeric(survT < censorT)
time    <- pmin(survT,censorT)

# Save the data and ground-truth coefficients to an .RData file.
save(list = c("b","geno","time","status"),
     file = "survival_demo.RData")
resaveRdaFiles("survival_demo.RData")
