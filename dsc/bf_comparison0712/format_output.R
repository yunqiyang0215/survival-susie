library(dscrutils)
res <- dscquery("dsc_result", c("simulate.censor_lvl", "simulate.allele_freq", 
                                "simulate.n", "simulate.b1", "compute_bf.result"))
saveRDS(res, "res_bf.rds")
