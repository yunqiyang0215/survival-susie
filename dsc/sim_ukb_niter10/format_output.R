# This file needs large memory to run....

library(dscrutils)

susie <- dscquery("dsc_result",
                  c("simulate.cor_type", "simulate.num_effect",  "simulate.censor_lvl",
                    "simulate.b", "simulate.is_effect", "susie.effect_estimate", "susie.pip", "susie.cs"), ignore.missing.files = TRUE)
saveRDS(susie, "susie.rds")

survsvb <- dscquery("dsc_result",
                    c("simulate.cor_type", "simulate.num_effect", "simulate.censor_lvl",
                      "simulate.b", "simulate.is_effect", "svb.effect_estimate", "svb.pip"))
saveRDS(survsvb, "survsvb.rds")


bvsnlp <- dscquery("dsc_result",
                   c("simulate.cor_type", "simulate.num_effect", "simulate.censor_lvl",
                     "simulate.b", "simulate.is_effect", "bvsnlp.effect_estimate", "bvsnlp.pip"))
saveRDS(bvsnlp, "bvsnlp.rds")


r2b <- dscquery("dsc_result",
                c("simulate.cor_type", "simulate.num_effect", "simulate.censor_lvl",
                  "simulate.b", "simulate.is_effect", "r2b.effect_estimate", "r2b.pip"), ignore.missing.files = TRUE)
saveRDS(r2b, "r2b.rds")


susie_rss <- dscquery("dsc_result",
                c("simulate.cor_type", "simulate.num_effect", "simulate.censor_lvl",
                  "simulate.b", "simulate.is_effect", "susie_rss.effect_estimate", "susie_rss.pip", "susie_rss.cs"))
saveRDS(susie_rss, "susie_rss.rds")


### Output time
time <- dscquery("dsc_result", c("susie.DSC_TIME", "svb.DSC_TIME", "bvsnlp.DSC_TIME",
                                 "susie_rss.DSC_TIME", "r2b.DSC_TIME"), ignore.missing.files = TRUE)
saveRDS(time, "time.rds")
