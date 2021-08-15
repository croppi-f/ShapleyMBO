###################################################### 
#########       getShapleyRes    #####################
######################################################
#- used to extract the results data frame fro the Shapley object
getShapleyRes = function(shapley, type) {
  # extract the usual iml::Shapley$results
  df = shapley$results %>% dplyr::mutate(
    pred.interest = shapley$y.hat.interest,
    pred.average = shapley$y.hat.average
   )
  
  # for the CB we need to extract some additional information to construct phi and phi.var later on.
  # - dataDesign: "the Design matrix after intervention" (see #61 on github repo iml/R/InterpretationMethod.R)
  # - qResults: "the quantity of interest from the black box model prediction" (see #65 on github repo iml/R/InterpretationMethod.R)
  # in case remove names of qResults
  res = switch(type,
               "simple" = df,
               "detailed" = list(df, 
                                 colnames(shapley$.__enclos_env__$private$dataDesign), 
                                 shapley$.__enclos_env__$private$qResults
                            )
               )
  return(res)
}

###################################################### 
#########       mergeShapleyRes     ##################
######################################################
#- used to merge the results of Shapley Mean and Shapley Se and to compute the
#  actual and average cb as well as the cb contributions with computePhiCb
#- this function actually genertaye the result data frame in each iteration
mergeShapleyRes = function(shapley.mean, shapley.se, lambda, max.mult, sample.size.s) {
  # 1. extract Shapely results 
  res.mean = lapply(shapley.mean, function(x) x[[1]]) %>% dplyr::bind_rows(.id = "iter")
  res.se = lapply(shapley.se, function(x) x[[1]]) %>% dplyr::bind_rows(.id = "iter")
  
  # 2. bind the results
  res = dplyr::bind_rows(res.mean, res.se, .id = "contribution") %>%
    dplyr::mutate(contribution = factor(contribution, labels=c("mean", "se")))
  
  # 3. reshape the results to create cb.interest and cb.average and scaled M and SE contributions
  res.wide = tidyr::pivot_wider(res, 
                                names_from = contribution, 
                                values_from = c(pred.interest, pred.average, phi, phi.var)) %>%
    dplyr::mutate(
      phi_se_scaled = -1 * lambda * phi_se, #scaling SE contributions with -1 * lambda
      phi_mean_scaled = max.mult * phi_mean, # scaling Mean contributions (needed for max problems, when max.mult = -1)
      pred.interest_cb = max.mult * pred.interest_mean - lambda * pred.interest_se,
      pred.average_cb = max.mult * pred.average_mean - lambda * pred.average_se,
      #phi_cb = max.mult * phi_mean - lambda * phi_se, # old version.  now, computed with computePhiCbLinear
    )
  names(res.wide$pred.average_cb) = rep(".prediction", nrow(res.mean))
  
  # 4. extract qResults to compute phi_cb and phi.var_cb with computePhiCb
  qR.mean = lapply(shapley.mean, function(x) x[[3]])
  qR.se = lapply(shapley.se, function(x) x[[3]])
  names.mean = lapply(shapley.mean, function(x) x[[2]]) %>% unique() %>% unlist()
  names.se = lapply(shapley.se, function(x) x[[2]]) %>% unique() %>% unlist()
  if(all.equal(names.mean, names.se)) names.dD = names.mean
  
  # 5. compute the SV of the CB using computePhiCb
  res.cb = mapply(computePhiCb, qR.mean, qR.se,
    MoreArgs = list(
      names = names.dD,
      lambda = lambda, max.mult = max.mult, sample.size.s = sample.size.s
    ),
    SIMPLIFY = FALSE, USE.NAMES = TRUE
  ) %>%
    dplyr::bind_rows(.id = "iter")
  
  # 6. bind res.wide and res.cb together
  res = dplyr::left_join(res.wide, res.cb, by = c("iter", "feature"))
  res = as.data.frame(res)
  
  return(res)
}


###################################################### 
#########       computePhiCb        ##################
######################################################
#- subfunction of mergeShapleyRes
#- this is the heart of the Shapley decomposition, where the CB contributions are
#  reconstructed using the Lineatiy Axiom
#- used to compute the Shapley Values of the LCB (phi_cb) and the variance of the estimation (phi.var_cb)
#- the main reason for this function is to correctly phi.var_cb 
#- Unlike cb. interest and average. to  compute phi_cb and phi.var_cb we need some 
#  additional private methods of the iml::Shapley object (see getShapleyRes type = "detailed")
#- for more details refer to the private method aggregate() in the iml repository R/Shapley.R on Github at line # 119 
#- computePhiCb is actually an extension of the method aggregate()
computePhiCb = function(qR.mean, qR.se, names, lambda, max.mult, sample.size.s) {
  
  # mean
  pred.with.k.mean = qR.mean[
    1:(nrow(qR.mean) / 2), , drop = FALSE
  ]
  pred.without.k.mean = qR.mean[
    (nrow(qR.mean) / 2 + 1) : nrow(qR.mean), , drop = FALSE
  ]
  
  # se
  pred.with.k.se = qR.se[
    1:(nrow(qR.se) / 2), , drop = FALSE
  ]
  pred.without.k.se = qR.se[
    (nrow(qR.se) / 2 + 1) : nrow(qR.se), , drop = FALSE
  ]
  
  # computing the cb
  pred.with.k = max.mult * pred.with.k.mean - lambda * pred.with.k.se
  pred.without.k =  max.mult * pred.without.k.mean - lambda * pred.without.k.se
  pred.diff = pred.with.k - pred.without.k
  
  cnames = colnames(pred.diff)
  
  pred.diff = cbind(
    data.table::data.table(feature = rep(names, times = sample.size.s)),
    pred.diff
  )
  
  pred.diff = data.table::melt(pred.diff, variable.name = "class", value.name = "value", measure.vars = cnames)
  res = pred.diff[, list("phi_cb" = mean(value), "phi.var_cb" = var(value)), by = c("feature", "class")] %>%
    dplyr::select(feature, phi_cb, phi.var_cb)
  
  return(res)
}
