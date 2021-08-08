source("R/plotShapleyCT.R")
source("R/plotShapleyCF.R")
###################################################### 
#########       plotShapleyMBO        ################
###################################################### 
#- two interesting view to see how the SV evolves (use argument by):
#  1) over the iterations (by = "time"): we do not expect any structure
#  2) with increasing AF values (by = "af"): we expect some trends
#- Note that in general AF decreases with increasing iters (e.g. LCB -> 0), so caution when using by= "af"
#- shapley.mbo the data frame containing the results of ShapleyAfMBO
#- iter.interest is used to plot results of single iterations. If NULL all stored iterations all plotted
#  for BAR PLOTS only one iter can be shown
#- by is used for LINE PLOTS and XPLXPL to display results with increasing iter of increasing
#  "prediction"
#- contribution = TRUE is used to call plotShapleyCBLinear(), which decompose the results in 
#  different contribution plots
#- decomp is used for plotShapleyCBLinear() only. Different decomposition plots are displayed depending
#  on the type argument. Possible values are 1,2,3.
#- scale is used by plotShapleyCbLinear.imp only. If = TRUE importance scores are scaled between 0 and 1.
#  Otherwise remain unscaled.

# TO DO: 
# - round the feature values to 4 digits when plotting
#%>% tidyr::separate(feature.value, into = c("f", "feature.value"), sep ="=", convert = TRUE)
# # rounding the feature value
# df$feature.value = round(df$feature.value, 4)
# df = df %>% tidyr::unite(f, feature.value, col = "feature.value", sep = "=", remove = TRUE)
# df = as.data.frame(df)

plotShapleyMBO = function(shapley.mbo, infill.mbo, lambda.mbo, sample.size = NULL, type, iter.interest = NULL,
                          avg.phi = FALSE, ci = FALSE, ci.alpha = 0.05, decomp = 1L) {
  # assertions & checks
  checkmate::assertDataFrame(shapley.mbo, any.missing = FALSE, all.missing = FALSE)
  checkmate::assertChoice(infill.mbo, c("mean", "se", "ei", "cb", "eqi", "aei", "adacb"))
  checkmate::assertNumber(lambda.mbo, null.ok = TRUE)
  checkmate::assertChoice(type, c("line", "bar"))
  checkmate::assertLogical(avg.phi, any.missing = FALSE, len = 1)
  checkmate::assertLogical(ci, any.missing = FALSE, len = 1)
  checkmate::assertChoice(ci.alpha, c(0.01, 0.05, 0.1))
  checkmate::assertNumber(sample.size, lower = 0, null.ok = TRUE)
  checkmate::assertChoice(decomp, c(1L, 2L, 3L))
  
  if(type == "bar" && (is.null(iter.interest) || length(iter.interest) > 1))
    stop("Bar-plot only available for one iteration at once")
  
  if(type == "line" && length(iter.interest) == 1)
    stop("Line Plot for single iteration not useful, use type = 'bp' instead")
  
  if(avg.phi == TRUE && ci == TRUE)
    stop("Choose between avg.phi or ci, but not both")
  
  if(ci == TRUE && is.null(sample.size))
    stop("Provide sample size used in ShapleyMBO")
  
  if(ci == TRUE && is.null(lambda.mbo))
    stop("Provide cb lambda used in BO run")
  
  if(avg.phi == TRUE && type != "line")
    stop("Average phi only for line plots available")

  #make numeric for line plots
  shapley.mbo$iter = as.numeric(shapley.mbo$iter)
  # only stored iterations can be plotted
  if(is.null(iter.interest)) {
    stored = sort(unique(shapley.mbo$iter))
  } else{
    iter.interest = as.numeric(iter.interest)
    checkmate::assertSubset(iter.interest, unique(shapley.mbo$iter))
    stored = sort(unique(iter.interest))
    if(length(stored) > 1 && max(stored) != length(stored)) {
      warning("Not adiacent iters have been selected. Line plots could be misleading.")
    }
  }
  # select the interesting iters
  shapley.mbo = shapley.mbo[shapley.mbo$iter %in% stored, ]
  
  if("phi_mean" %in% colnames(shapley.mbo)) {
    plot = plotShapleyCT(data = shapley.mbo, type = type, avg.phi = avg.phi, ci = ci, ci.alpha = ci.alpha, sample.size = sample.size, lambda = lambda.mbo, decomp = decomp)
  } else {               
    plot = plotShapleyCF(data = shapley.mbo, infill.mbo = infill.mbo, type = type, avg.phi = avg.phi, ci = ci, ci.alpha = ci.alpha, sample.size = sample.size)
  }                      
  
  return(plot)
}

