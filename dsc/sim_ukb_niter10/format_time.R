library(dscrutils)

susie <- dscquery("dsc_result",
                  c("simulate.cor_type", "simulate.num_effect", "simulate.censor_lvl", "susie.DSC_TIME")
                  ,ignore.missing.files = TRUE)

bvsnlp <- dscquery("dsc_result",
                   c("simulate.cor_type", "simulate.num_effect", "simulate.censor_lvl", "bvsnlp.DSC_TIME")
                   ,ignore.missing.files = TRUE)


susie_rss <- dscquery("dsc_result",
                      c("simulate.cor_type", "simulate.num_effect", "simulate.censor_lvl", "susie_rss.DSC_TIME")
                      ,ignore.missing.files = TRUE)

svb <- dscquery("dsc_result",
                c("simulate.cor_type", "simulate.num_effect", "simulate.censor_lvl", "svb.DSC_TIME")
                ,ignore.missing.files = TRUE)

r2b <- dscquery("dsc_result",
                c("simulate.cor_type", "simulate.num_effect", "simulate.censor_lvl", "r2b.DSC_TIME")
                ,ignore.missing.files = TRUE)


res = list(susie, bvsnlp, susie_rss, svb, r2b)
saveRDS(res, "time_comparison.rds")

