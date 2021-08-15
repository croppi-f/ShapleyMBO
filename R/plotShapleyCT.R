library(patchwork)
###################################################### 
#########       plotShapleyCT        #################
###################################################### 
# this is a subfunction of plotShapleyMBO to display the results
# of ShapleyMBO when contribution = TRUE.

plotShapleyCT = function(data, type, avg.phi, ci, ci.alpha, sample.size, lambda, decomp) {
  plot = switch(type,
                "line" = plotShapleyCT.line(data = data, avg.phi = avg.phi, ci = ci, ci.alpha = ci.alpha, sample.size = sample.size, lambda = lambda, decomp = decomp),
                "bar" = plotShapleyCT.bar(data = data, ci = ci, ci.alpha = ci.alpha, sample.size = sample.size, lambda = lambda, decomp = decomp)
  )
  return(plot)
}
###################################################### 
#########   plotShapleyCT.bar         ################
###################################################### 
# Bar Plot in single iterations, equivalent to iml::plot(Shapley)
plotShapleyCT.bar = function(data, ci, ci.alpha, sample.size, lambda, decomp) {
  # for better visualization we round the feature value 
  data = data %>% tidyr::separate(feature.value, into = c("f", "feature.value"), sep ="=", convert = TRUE)
  data$feature.value = round(data$feature.value, 4)
  data = data %>% tidyr::unite(f, feature.value, col = "feature.value", sep = "=", remove = TRUE)
  data = as.data.frame(data)
  
  # determine the order of the parameters
  order = forcats::fct_reorder(data$feature.value, data$phi_cb)
  lev = levels(order)
  
  # theme of the plot
  theme = theme_bw() + 
    theme(text = element_text(size = 15),# size 15 for better axis reading
          #axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "bottom"
    )
  
  # reshaping the data to display the decomposition plot
  data.long = data %>% dplyr::select(iter, feature.value, mean = phi_mean, se = phi_se_scaled) %>%
    tidyr::pivot_longer(cols = c(mean, se), 
                        names_to = "contribution",
                        names_ptypes = list(contribution = factor()),
                        values_to = "phi")
  
  if(ci == TRUE) {
    # scaling variance of the (unscaled) se contributions
    data.long.ci = data %>% dplyr::select(iter, feature.value, mean = phi.var_mean, se = phi.var_se) %>%
      dplyr::mutate(se = lambda ^ 2 * se) %>% # scale the variance of phi_se (not sclaed) with lambda^2 
      tidyr::pivot_longer(cols = c(mean, se), 
                          names_to = "contribution",
                          names_ptypes = list(contribution = factor()),
                          values_to = "phi.var")
    #adding phi_var to the data.long set
    data.long.ci = dplyr::left_join(data.long, data.long.ci, by = c("iter", "feature.value", "contribution"))
    # computing the range of the CI
    data.long.ci$phi.var = sqrt(data.long.ci$phi.var / sample.size)
    data.long.ci$range = data.long.ci$phi.var * qt(1 - ci.alpha / 2, sample.size - 1) # we use t distr since var unknown
    
    #- since mean contr bar is stacked on top of se contr bar, adjust phi_mean in order
    # to make it compatible with sd range, otherwise y of geom_errorbar for phi_mean will not be centered correctly
    #- when m and se contributions have opposite signs the problem does not occur
    #- if they do, phi_mean should be as "heigh" as phi_cb
    data$phi_mod = ifelse(sign(data$phi_mean) == sign(data$phi_se_scaled), data$phi_cb, data$phi_mean)
    data.mod.long = dplyr::select(data, iter, feature.value, se = phi_se_scaled, mean = phi_mean, phi_mod) %>%
      tidyr::pivot_longer(cols = c(mean, se), 
                          names_to = "contribution",
                          names_ptypes = list(contribution = factor()),
                          values_to = "phi")
    # replace values of phi_mod where contribution = se with their original values
    # since phi_se_scaled is not affected by this problem (phi_se_scaled is the lower bar)
    data.mod.long[seq(2, nrow(data.mod.long), 2), "phi_mod"] = data.mod.long[seq(2, nrow(data.mod.long), 2), "phi"]
    
    #bind data.long.ci and data.mod.long together and get final dat set for the plot
    data.long.ci = dplyr::left_join(data.mod.long, data.long.ci, by = c("iter", "feature.value", "contribution", "phi"))
    
    # contribution plot
    plot.contr = ggplot(data = data.long.ci, aes(x = forcats::fct_relevel(feature.value, lev), fill = contribution)) +
      geom_col(aes(y = phi)) +
      geom_errorbar(aes(ymin =  phi_mod - range, ymax =  phi_mod + range), width = 0.5) + 
      geom_hline(yintercept = 0, alpha = 0.5, linetype = "longdash") +
      coord_flip() +
      labs(x = "", y = "phi (cb)", title = sprintf("cb decomposed (iter %s) \nActual: %.4f, Average: %.4f", data$iter, data$pred.interest_cb, data$pred.average_cb))
    
    if(decomp == 1) {plot = plot.contr + theme}
    if(decomp == 2) {
      # CB Plot
      # compute the range for the ci
      data$phi.var_cb = sqrt(data$phi.var_cb / sample.size)
      data$range.cb = data$phi.var_cb * qt(1 - ci.alpha / 2, sample.size - 1) # we use t distr since var unknown
      
      plot.cb = ggplot(data = data, aes(y = phi_cb, x = forcats::fct_relevel(feature.value, lev))) +
        geom_col(fill = "#999999") +
        geom_errorbar(aes(ymin = phi_cb - range.cb, ymax = phi_cb + range.cb), width = 0.5) + 
        geom_hline(yintercept = 0, alpha = 0.5, linetype = "longdash") +
        coord_flip() +
        theme +
        labs(x = "", y = "phi (cb)",
             title = sprintf("Actual cb: %.4f\nAverage cb: %.4f",data$pred.interest_cb, data$pred.average_cb))
      
      # define the limits of the plots s.t. all plots have the same scale
      all.phi = c(data$phi_cb - data$range.cb, 
                  data$phi_cb + data$range.cb, 
                  data.long.ci$phi_mod - data.long.ci$range, 
                  data.long.ci$phi_mod + data.long.ci$range
      )
      lim = c(min(all.phi), max(all.phi))
      
      plot = plot.contr / plot.cb &
        theme &
        scale_y_continuous(limits = lim)
    }
    if(decomp == 3) {
      # CB Plot
      # compute the range for the ci
      data$phi.var_cb = sqrt(data$phi.var_cb / sample.size)
      data$range.cb = data$phi.var_cb * qt(1 - ci.alpha / 2, sample.size - 1) # we use t distr since var unknown
      plot.cb = ggplot(data = data, aes(y = phi_cb, x = forcats::fct_relevel(feature.value, lev))) +
        geom_col(fill = "#999999") +
        geom_errorbar(aes(ymin = phi_cb - range.cb, ymax = phi_cb + range.cb), width = 0.5) + 
        geom_hline(yintercept = 0, alpha = 0.5, linetype = "longdash") +
        coord_flip() +
        theme +
        labs(x = "", y = "phi (cb)",
             title = sprintf("Actual cb: %.4f\nAverage cb: %.4f",data$pred.interest_cb, data$pred.average_cb))
      # MEAN Plot
      # compute the range for the ci
      data$phi.var_mean = sqrt(data$phi.var_mean / sample.size)
      data$range.mean = data$phi.var_mean * qt(1 - ci.alpha / 2, sample.size - 1) # we use t distr since var unknown
      plot.mean = ggplot(data = data, aes(y = phi_mean, x = forcats::fct_relevel(feature.value, lev))) +
        geom_col(fill = "#999999") +
        geom_errorbar(aes(ymin = phi_mean - range.mean, ymax = phi_mean + range.mean), width = 0.5) + 
        geom_hline(yintercept = 0, alpha = 0.5, linetype = "longdash") +
        coord_flip() +
        theme +
        labs(x = "", y = "phi (m)",
             title = sprintf("Actual mean: %.4f\nAverage mean: %.4f",data$pred.interest_mean, data$pred.average_mean))
      
      # SE Plot (not scaled with lambda)
      # compute the range for the ci
      data$phi.var_se = sqrt(data$phi.var_se / sample.size) # here not scaled with lambda^2
      data$range.se = data$phi.var_se * qt(1 - ci.alpha / 2, sample.size - 1) # we use t distr since var unknown
      plot.se = ggplot(data = data, aes(y = phi_se, x = forcats::fct_relevel(feature.value, lev))) +
        geom_col(fill = "#999999") +
        geom_errorbar(aes(ymin = phi_se - range.se, ymax = phi_se + range.se), width = 0.5) + 
        geom_hline(yintercept = 0, alpha = 0.5, linetype = "longdash") +
        coord_flip() +
        theme +
        labs(x = "", y = "phi (se)",
             title = sprintf("Actual se: %.4f\nAverage se: %.4f",data$pred.interest_se, data$pred.average_se))
      
      # define the limits of the plots s.t. all plots have the same scale
      all.phi = c(data$phi_cb - data$range.cb, 
                  data$phi_cb + data$range.cb, 
                  data.long.ci$phi_mod - data.long.ci$range, 
                  data.long.ci$phi_mod + data.long.ci$range,
                  data$phi_mean - data$range.mean, 
                  data$phi_mean + data$range.mean,
                  data$phi_se - data$range.se, 
                  data$phi_se + data$range.se
      )
      lim = c(min(all.phi), max(all.phi))
      plot = plot.contr + plot.cb + plot.se + plot.mean +
        plot_layout(nrow = 2) &
        theme & 
        scale_y_continuous(limits = lim) 
    }
    
  }
  
  if(ci == FALSE) {
    # contribution plot
    plot.contr = ggplot(data = data.long, aes(x = forcats::fct_relevel(feature.value, lev), y = phi, fill = contribution)) +
      geom_col() +
      geom_hline(yintercept = 0, alpha = 0.5, linetype = "longdash") +
      coord_flip() +
      labs(x = "", y = "phi (cb)", title = sprintf("cb decomposed (iter %s) \nActual: %.4f, Average: %.4f", data$iter, data$pred.interest_cb, data$pred.average_cb))
    
    if(decomp == 1) {plot = plot.contr + theme}

    if(decomp == 2) {
      # CB Plot
      plot.cb = ggplot(data = data, aes(y = phi_cb, x = forcats::fct_relevel(feature.value, lev)))  +
        geom_col(fill = "#999999") +
        geom_hline(yintercept = 0, alpha = 0.5, linetype = "longdash") +
        coord_flip() +
        labs(x = "", y = "phi (cb)",
             title = sprintf("Actual cb: %.4f\nAverage cb: %.4f",data$pred.interest_cb, data$pred.average_cb))
      # define the limits of the plots s.t. all plots have the same scale
      
      all.phi = c(data$phi_cb, data$phi_mean, data$phi_se_scaled,
                  data$phi_mean + data$phi_se_scaled)
      lim = c(min(all.phi), max(all.phi))

      plot = plot.contr / plot.cb &
        theme &
        scale_y_continuous(limits = lim)
    }

    if(decomp == 3) {
      # CB plot
      plot.cb = ggplot(data = data, aes(y = phi_cb, x = forcats::fct_relevel(feature.value, lev)))  +
        geom_col(fill = "#999999") +
        geom_hline(yintercept = 0, alpha = 0.5, linetype = "longdash") +
        coord_flip() +
        labs(x = "", y = "phi (cb)",
             title = sprintf("Actual cb: %.4f\nAverage cb: %.4f",data$pred.interest_cb, data$pred.average_cb))
      #Mean
      plot.mean = ggplot(data = data, aes(y = phi_mean, x = forcats::fct_relevel(feature.value, lev))) +
        geom_col(fill = "#999999") +
        geom_hline(yintercept = 0, alpha = 0.5, linetype = "longdash") +
        coord_flip() +
        labs(x = "", y = "phi (m)",
             title = sprintf("Actual mean: %.4f\nAverage mean: %.4f", data$pred.interest_mean, data$pred.average_mean)
        )
      #Se
      plot.se = ggplot(data = data, aes(y = phi_se, x = forcats::fct_relevel(feature.value, lev))) +
        geom_col(fill = "#999999") +
        geom_hline(yintercept = 0, alpha = 0.5, linetype = "longdash") +
        coord_flip() +
        xlab("") +
        labs(x = "", y = "phi (se)",
             title = sprintf("Actual se: %.4f \nAverage se: %.4f", data$pred.interest_se, data$pred.average_se)
        )
      # define the limits of the plots s.t. all plots have the same scale
      all.phi = c(data$phi_mean, data$phi_se, data$phi_cb, data$phi_se_scaled, 
                  data$phi_mean + data$phi_se_scaled)
      lim = c(min(all.phi), max(all.phi))
      plot = plot.contr + plot.cb + plot.se + plot.mean +
        plot_layout(nrow = 2) &
        theme & 
        scale_y_continuous(limits = lim) 
    }
  }
  
  return(plot)
}

###################################################### 
#########   plotShapleyCT.line         ###############
###################################################### 
plotShapleyCT.line = function(data, avg.phi, ci, ci.alpha, sample.size, lambda, decomp) {
  # determine the order of the parameters
  order = forcats::fct_reorder(data$feature, data$phi_cb)
  
  # theme of the plot
  theme = theme_bw() + 
    theme(text = element_text(size = 15),
          legend.position = "bottom"
    )
  
  if(avg.phi == TRUE) {
    # compute the average phi for all the contributions
    avg.phi = data %>% dplyr::select(iter, feature, phi_mean, phi_se_scaled, phi_se, phi_cb) %>%
      dplyr::group_by(feature) %>% 
      dplyr::summarise(avg.mean = mean(phi_mean), avg.se_scaled = mean(phi_se_scaled), avg.se = mean(phi_se), avg.cb = mean(phi_cb))
    data = dplyr::left_join(data, avg.phi, by = "feature")
    
    # Cb Plot
    plot.cb = ggplot(data = data, aes(x = iter, group = feature)) +
      geom_hline(yintercept = 0, alpha = 0.5, linetype = "longdash") +
      geom_line(aes(y = phi_cb), color =  "#999999", alpha = 0.5) +
      geom_line(aes(y = avg.cb, color = order)) +
      theme +
      labs(x = "iter", y = "phi (cb)", color = "feature")
    # Mean Plot
    plot.mean = ggplot(data = data, aes(x = iter, group = feature)) +
      geom_hline(yintercept = 0, alpha = 0.5, linetype = "longdash") +
      geom_line(aes(y = phi_mean), color =  "#999999", alpha = 0.5) +
      geom_line(aes(y = avg.mean, color = order)) +
      theme + 
      labs(x = "iter", y = "phi (m)", color = "feature")
    # Se Plot
    plot.se = ggplot(data = data, aes(x = iter, group = feature)) +
      geom_hline(yintercept = 0, alpha = 0.5, linetype = "longdash") +
      geom_line(aes(y = phi_se_scaled), color =  "#999999", alpha = 0.5) +
      geom_line(aes(y = avg.se_scaled, color = order)) +
      theme + 
      labs(x = "iter", y = "phi (se)", color = "feature")
    
    if(decomp == 1) {
      # same scale for the plots
      all.phi = c(data$phi_cb,
                  data$phi_mean, 
                  data$phi_se_scaled
      )
      lim = c(min(all.phi), max(all.phi))
      plot = plot.cb + (plot.mean / plot.se) +
        plot_layout(guides = "collect") & 
        theme(legend.position = "bottom")  & 
        scale_y_continuous(limits = lim)
    }
    
    if(decomp == 2) {
      plot = plot.mean
    }
    
    if(decomp == 3) {
      # Se Plot (not scaled)
      plot.se.not = ggplot(data = data, aes(x = iter, group = feature)) +
        geom_hline(yintercept = 0, alpha = 0.5, linetype = "longdash") +
        geom_line(aes(y = phi_se), color =  "#999999", alpha = 0.5) +
        geom_line(aes(y = avg.se, color = order)) +
        theme + 
        labs(x = "iter", y = "phi (se) not scaled", color = "feature")
      
      all.phi = c(data$phi_se_scaled,
                  data$phi_se
      )
      lim = c(min(all.phi), max(all.phi))
      
      plot = plot.se + plot.se.not +
      plot_layout(guides = "collect") & 
        theme(legend.position = "bottom")  & 
        scale_y_continuous(limits = lim)
    }
  }
  
  if(ci == TRUE) {
    # compute the range for the CI for contributions cb, mean and se
    data$phi.var_mean = sqrt(data$phi.var_mean / sample.size)
    data$range.mean = data$phi.var_mean * qt(1 - ci.alpha / 2, sample.size - 1) # we use t distr since var unknown
    data$phi.var_se_scaled = sqrt(lambda ^ 2 * (data$phi.var_se / sample.size)) # scale the variance with lambda^2 here
    data$range.se_scaled = data$phi.var_se_scaled * qt(1 - ci.alpha / 2, sample.size - 1)
    data$phi.var_cb = sqrt(data$phi.var_cb / sample.size)
    data$range.cb = data$phi.var_cb * qt(1 - ci.alpha / 2, sample.size - 1)
    # CB Plot
    plot.cb = ggplot(data = data, aes(x = iter, y = phi_cb, group = feature, color = order)) +
      geom_hline(yintercept = 0, alpha = 0.5, linetype = "longdash") +
      geom_line() +
      geom_errorbar(aes(ymin = phi_cb - range.cb, ymax = phi_cb + range.cb)) +
      theme + 
      labs(x = "iter", y = "phi (cb)", color = "feature")
    # Mean Plot
    plot.mean = ggplot(data = data, aes(x = iter, y = phi_mean, group = feature, color = order)) +
      geom_hline(yintercept = 0, alpha = 0.5, linetype = "longdash") +
      geom_line() +
      geom_errorbar(aes(ymin = phi_mean - range.mean, ymax = phi_mean + range.mean)) +
      theme + 
      labs(x = "iter", y = "phi (m)", color = "feature")
    # Se Scaled Plot
    plot.se = ggplot(data = data, aes(x = iter, y = phi_se_scaled, group = feature, color = order)) +
      geom_hline(yintercept = 0, alpha = 0.5, linetype = "longdash") +
      geom_line() +
      geom_errorbar(aes(ymin = phi_se_scaled - range.se_scaled, ymax = phi_se_scaled + range.se_scaled)) +
      theme + 
      labs(x = "iter", y = "phi (se)", color = "feature")
    
    if(decomp == 1) {
      #same scale for all plots
      all.phi = c(data$phi_cb - data$range.cb, 
                  data$phi_cb + data$range.cb, 
                  data$phi_mean - data$range.mean, 
                  data$phi_mean + data$range.mean,
                  data$phi_se_scaled - data$range.se_scaled, 
                  data$phi_se_scaled + data$range.se_scaled
      )
      lim = c(min(all.phi), max(all.phi))
      
      plot = plot.cb + (plot.mean / plot.se) +
      plot_layout(guides = "collect") & 
      theme(legend.position = "bottom")  & 
      scale_y_continuous(limits = lim)
    }
    
    if(decomp == 2) { 
      plot = plot.mean
    }
    
    if(decomp == 3) {
      data$phi.var_se = sqrt(data$phi.var_se / sample.size) # here not sclaed with lambda^2
      data$range.se = data$phi.var_se * qt(1 - ci.alpha / 2, sample.size - 1)
      # Se Plot (not scaled)
      plot.se.not = ggplot(data = data, aes(x = iter, y = phi_se, group = feature, color = order)) +
        geom_hline(yintercept = 0, alpha = 0.5, linetype = "longdash") +
        geom_line() +
        geom_errorbar(aes(ymin = phi_se - range.se, ymax = phi_se + range.se)) +
        theme + 
        labs(x = "iter", y = "phi (se) not scaled", color = "feature")
      
      all.phi = c(
        data$phi_se_scaled - data$range.se_scaled, 
        data$phi_se_scaled + data$range.se_scaled,
        data$phi_se - data$range.se, 
        data$phi_se + data$range.se
      )
      lim = c(min(all.phi), max(all.phi))
      
      plot = plot.se + plot.se.not +
        plot_layout(guides = "collect") & 
        theme(legend.position = "bottom")  & 
        scale_y_continuous(limits = lim)
    }
  }
  
  if(avg.phi == FALSE && ci == FALSE) {
    # CB Plot
    plot.cb = ggplot(data = data, aes(x = iter, y = phi_cb, group = feature, color = order)) +
      geom_hline(yintercept = 0, alpha = 0.5, linetype = "longdash") +
      geom_line() +
      theme +
      labs(x = "iter", y = "phi (cb)", color = "feature") 
    # Mean Plot
    plot.mean = ggplot(data = data, aes(x = iter, y = phi_mean, group = feature, color = order)) +
      geom_hline(yintercept = 0, alpha = 0.5, linetype = "longdash") +
      geom_line() +
      theme + 
      labs(x = "iter", y = "phi (m)", color = "feature")
    # Se Scaled Plot
    plot.se = ggplot(data = data, aes(x = iter, y = phi_se_scaled, group = feature, color = order)) +
      geom_hline(yintercept = 0, alpha = 0.5, linetype = "longdash") +
      geom_line() +
      theme + 
      labs(x = "iter", y = "phi (se)", color = "feature")
    
    if(decomp == 1) {
      all.phi = c(data$phi_cb,
                  data$phi_mean,
                  data$phi_se_scaled
      )
      lim = c(min(all.phi), max(all.phi))
      
      plot = plot.cb + (plot.mean / plot.se) +
        plot_layout(guides = "collect") & 
        theme(legend.position = "bottom")  & 
        scale_y_continuous(limits = lim)
    }
    
    if(decomp == 2) {
      plot = plot.mean
    }
    
    if(decomp == 3) {
      # Se Plot (not scaled)
      plot.se.not = ggplot(data = data, aes(x = iter, y = phi_se, group = feature, color = order)) +
        geom_hline(yintercept = 0, alpha = 0.5, linetype = "longdash") +
        geom_line() +
        theme + 
        labs(x = "iter", y = "phi (se) not scaled", color = "feature")
      
      all.phi = c(
        data$phi_se_scaled,
        data$phi_se
      )
      lim = c(min(all.phi), max(all.phi))
      
      plot = plot.se + plot.se.not +
      plot_layout(guides = "collect") & 
        theme(legend.position = "bottom")  & 
        scale_y_continuous(limits = lim)
    }
  }
  
  return(plot)
}