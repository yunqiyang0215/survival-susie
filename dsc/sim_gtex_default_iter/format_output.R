library(dscrutils)

susie <- dscquery("dsc_result",
                  c("simulate.cor_type", "simulate.num_effect",  "simulate.censor_lvl",
                    "simulate.b", "simulate.is_effect", "susie.effect_estimate", "susie.pip", "susie.cs"))

svb <- dscquery("dsc_result",
                c("simulate.cor_type", "simulate.num_effect", "simulate.censor_lvl",
                  "simulate.b", "simulate.is_effect", "svb.effect_estimate", "svb.pip"))


bvsnlp <- dscquery("dsc_result",
                   c("simulate.cor_type", "simulate.num_effect", "simulate.censor_lvl",
                     "simulate.b", "simulate.is_effect", "bvsnlp.effect_estimate", "bvsnlp.pip"))

rss <- dscquery("dsc_result",
                c("simulate.cor_type", "simulate.num_effect", "simulate.censor_lvl",
                  "simulate.b", "simulate.is_effect", "rss.effect_estimate", "rss.pip", "rss.cs"))

r2b <- dscquery("dsc_result",
                c("simulate.cor_type", "simulate.num_effect", "simulate.censor_lvl",
                  "simulate.b", "simulate.is_effect", "r2b.effect_estimate", "r2b.pip"))

saveRDS(susie, "susie.rds")
saveRDS(svb, "svb.rds")
saveRDS(bvsnlp, "bvsnlp.rds")
saveRDS(rss, "rss.rds")
saveRDS(r2b, "r2b.rds")




