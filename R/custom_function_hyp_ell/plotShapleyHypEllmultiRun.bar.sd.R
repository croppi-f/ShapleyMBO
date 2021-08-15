library(patchwork)
#- this function is the equivalent and less flexible version of plotShapleyCT.bar 
#  to analyze the average results over multiple BO runs
#- Differences: since we are dealing with multi runs and the Hyper Ellipsoid
#    1) we plot the average phi
#    2) for the unceratinty we use the sd instead of the CI
#    3) we display the average distance from the optimum value instead of the average parameter 
#       value on the vertical axis
#    4) no options to choose which plot to display (all 4 plots + uncertainty estimates)
#- this function is used as a faster and more clean way to analyze the results 
#  of the HyperEllipsoid
#-it should therefore not be used for other purposes & should not be used by others (not "exported")
#- for more details see plotShapleyCT.bar
#- shapley.mbo is a data frame containing the ShapleyMBO results for multiple BO runs
plotShapleyHypEllmultiRun.bar.sd = function(shapley.mbo, lambda.mbo, iter.interest) {
  # theme of the plot
  theme = theme_bw() + 
    theme(text = element_text(size = 15),
          legend.position = "bottom"
    )
  
  # get the data
  data = shapley.mbo %>% dplyr::select(
    lambda, iter, feature.value, 
    phi_mean, pred.interest_mean, pred.average_mean,
    phi_se, phi_se_scaled, pred.interest_se, pred.average_se,
    phi_cb, pred.interest_cb, pred.average_cb
  ) %>% dplyr::slice(which(lambda == as.character(lambda.mbo) & iter == as.character(iter.interest)))
  # split the feature value column and compute the distance to optimum value, 
  data = data %>% tidyr::separate(feature.value, into = c("feature", "value"), sep ="=", convert = TRUE)
  data$value = abs(data$value) # because optimum in at 0 --> abs(data$value - 0)
  
  # compute the average results (phi, prediction interest and average) over the runs
  avg.data = data %>% group_by(lambda, feature, iter) %>%
    dplyr::summarise(across(c(starts_with("phi"), starts_with("pred."),"value"),mean))
  # std.dev used as uncertainty of the estimate
  sd.phi = data %>% group_by(lambda, feature, iter) %>%
    dplyr::summarise(across(starts_with("phi"), sd))
  data = dplyr::left_join(avg.data, sd.phi, by = c("lambda", "feature", "iter"), suffix = c(".avg", ".sd"))
  
  # merge parameter and parameter value back together
  data$value = round(data$value, 4) # four digits for better visibility
  data = data %>% tidyr::unite(feature, value, col = "feature.value", sep = "=", remove = TRUE)
  data = as.data.frame(data)
  # determine the order of the features
  order = forcats::fct_reorder(data$feature.value, data$phi_cb.avg)
  lev = levels(order)
  
  # since mean contr bar is stacked on top of se contr bar, adjust phi_mean in order
  # to make it compatible with sd range
  data$phi_mod.avg = ifelse(
    sign(data$phi_mean.avg) == sign(data$phi_se_scaled.avg),
    data$phi_cb.avg, # phi_mean will be equal to phi_cb
    data$phi_mean.avg
  )
  
  # 1.Plot: contribution plot
  # transform to long format separately phi and sd
  data.phi.long = dplyr::select(data, lambda, feature.value, iter, mean = phi_mean.avg, se = phi_se_scaled.avg) %>% 
    tidyr::pivot_longer(cols = c(mean, se), names_to = "contribution", names_ptypes = list(contribution = factor()),
                        values_to = "phi")
  data.sd.long = dplyr::select(data, lambda, feature.value, iter, mean = phi_mean.sd, se = phi_se_scaled.sd) %>% 
    tidyr::pivot_longer(cols = c(mean, se), names_to = "contribution", names_ptypes = list(contribution = factor()),
                        values_to = "sd")
  # merge them
  data.sd.long = dplyr::left_join(data.phi.long, data.sd.long, by = c("lambda", "iter", "feature.value", "contribution"))
  
  # transform to long format data including phi_mod column
  data.mod.long = dplyr::select(data, iter, feature.value, se = phi_se_scaled.avg, mean = phi_mean.avg, phi_mod = phi_mod.avg) %>%
    tidyr::pivot_longer(cols = c(mean, se), names_to = "contribution", names_ptypes = list(contribution = factor()),
                        values_to = "phi")
  # replace values of phi_mod where contribution = se with their original values
  # since phi_se_scaled is not affected by this problem (phi_se_scaled is the lower bar)
  data.mod.long[seq(2, nrow(data.mod.long), 2), "phi_mod"] = data.mod.long[seq(2, nrow(data.mod.long), 2), "phi"]
  #join the data and obtain the final data set used for the plot
  data.contr.plot = dplyr::left_join(data.mod.long, data.sd.long, by = c("iter", "feature.value", "contribution", "phi"))
  # this is a subfunction used to customize the parameter axis labels, it replaces "x" with expression (theta)
  changeXinTheta = function(old.breaks) {
    split = strsplit(old.breaks, "=")
    new.breaks = list()
    for(i in seq_along(split)) {
      feature = split[[i]][1]
      value = split[[i]][2]
      if(feature == "x1") new.breaks[[i]] = expr(paste(theta[1], sep = "=", !!value))
      if(feature == "x2") new.breaks[[i]] = expr(paste(theta[2], sep = "=", !!value))
      if(feature == "x3") new.breaks[[i]] = expr(paste(theta[3], sep = "=", !!value))
      if(feature == "x4") new.breaks[[i]] = expr(paste(theta[4], sep = "=", !!value))
    }
    new.breaks
  }
  # contribution plot (decomposition plot)
  plot.contr = ggplot(data = data.contr.plot, aes(x = forcats::fct_relevel(feature.value, lev), fill = contribution)) +
    geom_col(aes(y = phi)) +
    geom_errorbar(aes(ymin =  phi_mod - sd, ymax =  phi_mod + sd), width = 0.5) + 
    geom_hline(yintercept = 0, alpha = 0.5, linetype = "longdash") +
    scale_x_discrete(labels = changeXinTheta)+
    coord_flip() +
    labs(x = "", y = "phi (cb)", title = sprintf("cb decomposed (iter %s) \nActual: %.4f, Average: %.4f", data$iter, data$pred.interest_cb, data$pred.average_cb))
  
  # 2.Plot: cb
  plot.cb = ggplot(data = data, aes(y = phi_cb.avg, x = forcats::fct_relevel(feature.value, lev))) +
    geom_col(fill = "#999999") +
    geom_errorbar(aes(ymin = phi_cb.avg - phi_cb.sd, ymax = phi_cb.avg + phi_cb.sd), width = 0.5) + 
    geom_hline(yintercept = 0, alpha = 0.5, linetype = "longdash") +
    scale_x_discrete(labels = changeXinTheta) +
    coord_flip() +
    theme +
    labs(x = "", y = "phi (cb)",
         title = sprintf("Actual cb: %.4f\nAverage cb: %.4f",data$pred.interest_cb, data$pred.average_cb))
  # 3.Plot: mean
  plot.mean = ggplot(data = data, aes(y = phi_mean.avg, x = forcats::fct_relevel(feature.value, lev))) +
    geom_col(fill = "#999999") +
    geom_errorbar(aes(ymin = phi_mean.avg - phi_mean.sd, ymax = phi_mean.avg + phi_mean.sd), width = 0.5) + 
    geom_hline(yintercept = 0, alpha = 0.5, linetype = "longdash") +
    scale_x_discrete(labels = changeXinTheta) +
    coord_flip() +
    theme +
    labs(x = "", y = "phi (mean)",
         title = sprintf("Actual mean: %.4f\nAverage mean: %.4f",data$pred.interest_mean, data$pred.average_mean))
  # 4.Plot: se
  plot.se = ggplot(data = data, aes(y = phi_se.avg, x = forcats::fct_relevel(feature.value, lev))) +
    geom_col(fill = "#999999") +
    geom_errorbar(aes(ymin = phi_se.avg - phi_se.sd, ymax = phi_se.avg + phi_se.sd), width = 0.5) + 
    geom_hline(yintercept = 0, alpha = 0.5, linetype = "longdash") +
    scale_x_discrete(labels = changeXinTheta) +
    coord_flip() +
    theme +
    labs(x = "", y = "phi (se)",
         title = sprintf("Actual se: %.4f\nAverage se: %.4f",data$pred.interest_se, data$pred.average_se))
  
  # define the limits of the plots s.t. all plots have the same range
  all.phi = c(data$phi_cb.avg - data$phi_cb.sd, 
              data$phi_cb.avg + data$phi_cb.sd, 
              data.contr.plot$phi_mod - data.contr.plot$sd, 
              data.contr.plot$phi_mod + data.contr.plot$sd,
              data$phi_mean.avg - data$phi_mean.sd, 
              data$phi_mean.avg+ data$phi_mean.sd,
              data$phi_se.avg - data$phi_se.sd, 
              data$phi_se.avg + data$phi_se.sd
  )
  lim = c(min(all.phi), max(all.phi))
  plot = plot.contr + plot.cb + plot.se + plot.mean +
    plot_layout(nrow = 2) &
    theme & 
    scale_y_continuous(limits = lim) 
  
  return(plot)
}
