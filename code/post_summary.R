calculate_tpr_vs_fdr <- function(pip, is_effect, ts){
  res <- matrix(NA, nrow = length(ts), ncol = 2)
  colnames(res) = c("tpr", "fdr")
  for (i in 1:length(ts)){
    pred_pos = pip >= ts[i]
    tp = pip >= ts[i] & is_effect == 1
    fp = pip >= ts[i] & is_effect == 0
    tpr = sum(tp)/sum(is_effect)
    fdr = sum(fp)/sum(pred_pos)
    res[i, ] = c(tpr, fdr)
  }
  return(res)
}



# coverage: the proportion of CSs that contain an effect variable
# @param dat_indx: the indx for the data from dsc
# @param res.cs: credible sets from dsc
calculate_cs_coverage = function(res.cs, res.is_effect, dat_indx){
  contain_status = c()
  for (indx in dat_indx){
    cs = res.cs[[indx]]$cs
    true_effect = which(res.is_effect[[indx]] >= 1)
    if (!is.null(cs)){
      for (j in 1:length(cs)){
        res = ifelse(sum(true_effect %in% unlist(cs[j])) ==  1, 1, 0)
        contain_status = c(contain_status, res)
      }
    }
  }
  coverage = sum(contain_status)/length(contain_status)
  return(coverage)
}

# @param res.cs: credible sets from dsc
# @param dat_indx: the indx for the data from dsc
# @p: number of variables in each simulation replicate.
get_cs_effect = function(res.cs, dat_indx, p){
  cs_effect = c()
  for (indx in dat_indx){
    effect = rep(0, p)
    cs_effect_indx = c(unlist(res.cs[[indx]]$cs))
    effect[cs_effect_indx] = 1
    cs_effect = c(cs_effect, effect)
  }
  return(cs_effect)
}
