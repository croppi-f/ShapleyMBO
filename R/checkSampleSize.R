################################################################################
#####################   checkSampleSize   ######################################
################################################################################
#'@description This function is used to check if the sample size used in ShapleyMBO is sufficiently
#'  high. For more details see Appendix C in the thesis paper.
#'@param shapley.mbo [\code{data.frame}] \cr
#'  \code{ShapleyMBO} results for one single iteraions only.
#'@return A string indicating if the sample for cb, mean and se contribution is high enough.

checkSampleSize = function(shapley.mbo) {
  # repeat the same steps for cb, mean and se: 
  # Cb
  # compute the payout
  payout.cb = unique(shapley.mbo$pred.interest_cb) - unique(shapley.mbo$pred.average_cb)
  # compute the sum of the sv
  sum.phi.cb = sum(shapley.mbo$phi_cb)
  # compute the efficiency error (gap)
  gap.cb = abs(sum.phi.cb - payout.cb)
  # compute the minimum distance between the phi (threshold)
  diffs.cb = outer(shapley.mbo$phi_cb, shapley.mbo$phi_cb, "-")
  threshold.cb = min(
    abs( 
      diffs.cb[upper.tri(diffs.cb)]
    )
  )
  # sample size check
  check.cb = ifelse(gap.cb < threshold.cb, "TRUE", "FALSE")
  
  #Mean
  payout.mean = unique(shapley.mbo$pred.interest_mean) - unique(shapley.mbo$pred.average_mean)
  sum.phi.mean = sum(shapley.mbo$phi_mean)
  gap.mean = abs(sum.phi.mean - payout.mean)
  diffs.mean = outer(shapley.mbo$phi_mean, shapley.mbo$phi_mean, "-")
  threshold.mean = min(
    abs(
      diffs.mean[upper.tri(diffs.mean)]
    )
  )
  check.mean = ifelse(gap.mean < threshold.mean, "TRUE", "FALSE")
  
  #Se
  payout.se = unique(shapley.mbo$pred.interest_se) - unique(shapley.mbo$pred.average_se)
  sum.phi.se = sum(shapley.mbo$phi_se)
  gap.se = abs(sum.phi.se - payout.se)
  diffs.se = outer(shapley.mbo$phi_se, shapley.mbo$phi_se, "-")
  threshold.se = min(
    abs(
      diffs.se[upper.tri(diffs.se)]
    )
  )
  check.se = ifelse(gap.se < threshold.se, "TRUE", "FALSE")
  
  return(sprintf("sample size high enough for: cb %s, mean %s, se %s", check.cb, check.mean, check.se))
}
