source("R/plotShapleyCT.R")
source("R/plotShapleyCF.R")
###################################################### 
#########       plotShapleyMBO        ################
###################################################### 
#'@description Displays the results of \code{ShapleyMBO}
#'
#'@details Choose to display single iterations results with bar
#'plot or sol called desirability paths. To arrange plots together the \pkg{patchwork}
#'is used. Individual plots can also be extracted from the stored object by selecting
#'the desired element as if they were list elements using [[]] brackets, e.g. [[1]] to get the first plot.
#'
#'@param shapley.mbo [\code{data.frame}] \cr
#'  ShapleyMBO result object
#'@param infill.mbo [\code{character(1)}] \cr
#'  Acquisition function used in the mbo optimization
#'@param lambda.mbo [\code{numeric(1)}] \cr
#'  Control parameter of the CB. Only necessary when the CB is used in the optimization.
#'  Default is \code{NULL}.
#'@param sample.size [\code{numeric(1)}] \cr
#'  The number of Monte Carlo samples used by ShapleyMBO. This argument is used to compute 
#'  the confidence intervals of the Shapley values.
#'@param type [\code{character(1)}] \cr
#'  Which type of plots should be generated? Choose \code{type = "bar"} for on individual iteration
#'  or \code{type = "line"} to plot desirability paths.
#'@param iter.interest [\code{numeric}] \cr
#'  Iteration for which the results should be plotted. Multiple are possible (numeric vector). 
#'  Default is \code{NULL}, which means that all possible iterations are plotted. Do not use \code{iter.interest = NULL}
#'  with \code{type = "bar"} because it generates an error. 
#'@param avg.phi [\code{logical(1)}] \cr
#'  Only possible in combination with \code{type = "line"}. This option plots the disability paths
#'  (grey curves) and highlights in the same plot the average contribution of parameters.
#'@param ci [\code{logical(1)}] \cr
#'  Should the confidence intervals of the contributions be plotted? Default is \code{FALSE}. Works for both 
#'  plot types and does not work if \code{avg.phi = TRUE}.
#'@param ci.alpha [\code{logical(1)}] \cr
#'  Significance level of the confidence interval. Three options available: 0.01, 0.05, 0.1. Default is \code{ci.alpha = 0.05}
#'@param decomp [\code{numeric(1)}] \cr
#'  Is used when \code{infill.mbo = "cb"} and allows to choose different plot arrangements depending on \code{type}. 
#'  For \code{type = "bar"}: \code{decomp = 1L} only decomposition plot, \code{decomp = 2L} decomposition and overall
#'  desraibility (cb contributions), \code{decomp = 3L} decomposition plot, cb, mean and se contribution.
#'  For \code{type = "line"}: \code{decomp = 1L} cb, mean and se contributions, \code{decomp = 2L} only
#'  mean contributions, \code{decomp = 3L} only se contributions, scaled and not scaled.
#'  
#'@return A plot object (list) of class \code{patchwork, gg, ggplot}.

plotShapleyMBO = function(shapley.mbo, infill.mbo, lambda.mbo, sample.size = NULL, type, iter.interest = NULL,
                          avg.phi = FALSE, ci = FALSE, ci.alpha = 0.05, decomp = 1L) {
  # input checking
  checkmate::assertDataFrame(shapley.mbo, any.missing = FALSE, all.missing = FALSE)
  checkmate::assertChoice(infill.mbo, c("mean", "se", "ei", "cb", "eqi", "aei", "adacb"))
  checkmate::assertNumber(lambda.mbo, null.ok = TRUE)
  checkmate::assertNumber(sample.size, lower = 0, null.ok = TRUE)
  checkmate::assertChoice(type, c("line", "bar"))
  checkmate::assertNumeric(iter.interest, lower = 1, any.missing = FALSE,
                           all.missing = FALSE, min.len = 1, unique = TRUE,
                           null.ok = TRUE)
  checkmate::assertLogical(avg.phi, any.missing = FALSE, len = 1)
  checkmate::assertLogical(ci, any.missing = FALSE, len = 1)
  checkmate::assertChoice(ci.alpha, c(0.01, 0.05, 0.1))
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
  # if not specified, only stored iterations can be plotted
  if(is.null(iter.interest)) {
    stored = sort(unique(shapley.mbo$iter))
  } else{
    iter.interest = as.numeric(iter.interest)
    #iter.interest must be in shapley.mbo
    checkmate::assertSubset(iter.interest, unique(shapley.mbo$iter))
    stored = sort(unique(iter.interest))
    if(length(stored) > 1 && max(stored) != length(stored)) {
      warning("Not adiacent iters have been selected. Line plots could be misleading.")
    }
  }
  # select the data out of shpaley.mbo
  shapley.mbo = shapley.mbo[shapley.mbo$iter %in% stored, ]
  
  if("phi_mean" %in% colnames(shapley.mbo)) {
    #ShapleyMBO was called with contribution = TRUE
    plot = plotShapleyCT(data = shapley.mbo, type = type, avg.phi = avg.phi, ci = ci, ci.alpha = ci.alpha, sample.size = sample.size, lambda = lambda.mbo, decomp = decomp)
  } else {
    #ShapleyMBO was called with contribution = FALSE
    plot = plotShapleyCF(data = shapley.mbo, infill.mbo = infill.mbo, type = type, avg.phi = avg.phi, ci = ci, ci.alpha = ci.alpha, sample.size = sample.size)
  }                      
  
  return(plot)
}

