pip_calibration <- function(pips, is_effect){
  ts = seq(0, 0.9, by = 0.1)
  res = matrix(NA, ncol = 3, nrow = length(ts))
  colnames(res) = c("range.start", "expected", "empirical")
  for (i in 1:length(ts)){
    lo = ts[i]
    up = lo + 0.1
    indx = which(pips >= lo & pips <= up)
    res[i, 1] = lo
    res[i, 2] = mean(pips[indx])
    res[i, 3] = sum(is_effect[indx])/ length(indx)
  }
  # when there is no pip in the certain range
  res2 = na.omit(res)
  return(res2)
}

plot_calibration <- function(dat.calibration, main){
  plot(dat.calibration[,2], dat.calibration[,3], xlim = c(0, 1), ylim = c(0, 1),
       main= main, xlab="Expected", ylab="Observed")
  abline(a = 0, b = 1, col = "red")
}

compare_pip <- function(pip1, pip2, is_effect){
  res = data.frame(cbind(pip1, pip2, is_effect))
  return(res)
}

plot_pip <- function(res.pip, labs, main){
  pch = 16
  cex = ifelse(res.pip$is_effect == 1, 0.75, 1.2)

  # Plotting the grey circles first
  grey_indices <- res.pip$is_effect != 1
  plot(res.pip$pip1[grey_indices], res.pip$pip2[grey_indices], col = "dark grey", pch = 1, cex = cex[grey_indices], xlab = labs[1], ylab = labs[2], xlim = c(0, 1), ylim = c(0, 1), main = main)

  red_indices <- res.pip$is_effect == 1
  points(res.pip$pip1[red_indices], res.pip$pip2[red_indices], col = "red", pch = pch, cex = cex[red_indices])
}
