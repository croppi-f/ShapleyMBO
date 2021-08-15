######################################################
#########       plotShapleyCF        #################
######################################################
# this is a subfunction of plotShapleyMBO to display the results
# of ShapleyMBO when contribution = FALSE.

plotShapleyCF = function(data, infill.mbo, type, avg.phi, ci, ci.alpha, sample.size) {
  plot = switch(type,
                "line" = plotShapleyCF.line(data = data, infill.mbo = infill.mbo, avg.phi = avg.phi, ci = ci, ci.alpha = ci.alpha, sample.size = sample.size),
                "bar" = plotShapleyCF.bar(data = data, infill.mbo = infill.mbo, ci = ci, ci.alpha = ci.alpha, sample.size = sample.size),
  )
  return(plot)
}

######################################################
#########       plotShapleyCF.line   #################
######################################################
plotShapleyCF.line = function(data, infill.mbo, avg.phi, ci, ci.alpha, sample.size) {
  # theme of the plot
  theme = theme_bw() + 
    theme(text = element_text(size = 15),# size 15 for better axis reading
          #axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "bottom"
    )
  
  if(avg.phi == TRUE) {
    #compute the avg contribution
    avg.phi = data %>%
      dplyr::group_by(feature) %>%
      dplyr::summarise(avg.phi = mean(phi))
    data = dplyr::left_join(data, avg.phi, by = "feature")
    
    plot = ggplot(data = data, aes(x = iter, group = feature)) + 
      geom_hline(yintercept = 0, alpha = 0.5, linetype = "longdash") +
      geom_line(aes(y = phi), color =  "#999999", alpha = 0.5) +
      geom_line(aes(y = avg.phi, color = feature)) +
      theme + 
      labs(
        x = "iter", y = sprintf("phi (%s)", infill.mbo)
      )
  }
  if(ci == TRUE) {
    #compute the range of the CI:
    data$phi.var = sqrt(data$phi.var / sample.size)
    data$range = data$phi.var * qt(1 - ci.alpha / 2, sample.size - 1) # we use t distr since var unknown
    
    plot = ggplot(data = data, aes(x = iter, y = phi, group = feature, color = feature)) + 
      geom_hline(yintercept = 0, alpha = 0.5, linetype = "longdash") +
      geom_line() +
      geom_point() +
      geom_errorbar(aes(ymin = phi - range, ymax = phi + range), width = 0.5) +
      theme + 
      labs(
        x = "iter", y = sprintf("phi (%s)", infill.mbo)
      )
  }
  if(ci == FALSE && avg.phi == FALSE) {
    plot = ggplot(data = data, aes(x = iter, y = phi, group = feature, color = feature)) +
      geom_hline(yintercept = 0, alpha = 0.5, linetype = "longdash") +
      geom_line() +
      theme + 
      labs(
        x = "iter", y = sprintf("phi (%s)", infill.mbo)
      )
  }

  return(plot)
}

######################################################
#########       plotShapleyCF.bar    #################
######################################################
# Bar Plot in single iterations, equivalent to iml::plot(Shapley)
plotShapleyCF.bar = function(data, infill.mbo, ci, ci.alpha, sample.size) {
  # for better visualization we round the feature value to 4 digits
  data = data %>% tidyr::separate(feature.value, into = c("f", "feature.value"), sep ="=", convert = TRUE)
  data$feature.value = round(data$feature.value, 4)
  data = data %>% tidyr::unite(f, feature.value, col = "feature.value", sep = "=", remove = TRUE)
  data = as.data.frame(data)
  
  # defining the order of the parameters
  order = forcats::fct_reorder(data$feature.value, data$phi)
  lev = levels(order)
  
  # theme of the plot
  theme = theme_bw() + 
    theme(text = element_text(size = 15),# size 15 for better axis reading
          #axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "bottom"
    )
  
  if(ci == TRUE) {
    data$phi.var = sqrt(data$phi.var / sample.size)
    data$range = data$phi.var * qt(1 - ci.alpha / 2, sample.size - 1) # we use t distr since var unknown
    
    plot = ggplot(data = data, aes(y = phi, x = forcats::fct_relevel(feature.value, lev))) +
      geom_col(fill = "#999999") +
      geom_errorbar(aes(ymin = phi - range, ymax = phi + range), width = 0.5) + 
      geom_hline(yintercept = 0, alpha = 0.5, linetype = "longdash") +
      coord_flip() +
      theme +
      labs(x = "", y = sprintf("phi (%s)", infill.mbo), 
           title = sprintf(
             "iter %s \nActual %s: %.3f, Avg %s: %.3f",
             data$iter, infill.mbo, data$pred.interest, infill.mbo, data$pred.average
           )
      )
  }
  if(ci == FALSE) {
    plot = ggplot(data = data, aes(y = phi, x = forcats::fct_relevel(feature.value, lev))) +
      geom_col(fill = "#999999") +
      geom_hline(yintercept = 0, alpha = 0.5, linetype = "longdash") +
      coord_flip() +
      theme +
      labs(x = "", y = sprintf("phi (%s)", infill.mbo), 
        title = sprintf(
          "iter %s \nActual %s: %.3f, Avg %s: %.3f",
          data$iter, infill.mbo, data$pred.interest, infill.mbo, data$pred.average
        )
      )
  }
  
  return(plot)
}
