source("R/utils_ShapleyMBO.R")
source("R/ShapleyAf.R")
source("R/ShapleyMean.R")
source("R/ShapleySe.R")
source("R/PredictorAf.R")
source("R/utils_PredictorAf.R")
###################################################### 
#######     ShapleyMBO_mclapply     ##################
######################################################
#- this is an alternative version of ShapleyMBO which allows
#  parallel computation. Works only on MacOS and Linux, not on Windows
#- for a description of the function see "ShapleyMBO.R"
#- no.cores: Number of cores to use for the parallel computation. Default is 1.
ShapleyMBO_mclapply = function(res.mbo,
                               iter.interest = NULL, # explain instance proposed in iter.interest
                               sample.size = 100,
                               contribution = FALSE,
                               seed = 1,
                               no.cores = 1) {
  # Input checking before extracting infos from MBO object
  checkmate::assertClass(res.mbo, classes = c("MBOSingleObjResult", "MBOResult"))
  
  # Extract useful information from MBO object
  opdf =  as.data.frame(res.mbo$opt.path)
  #ParSet and Control object
  ps.mbo = res.mbo$opt.path$par.set
  ctrl.mbo = res.mbo$control
  # Pred Models
  sm = res.mbo$models
  af = ctrl.mbo$infill.crit$fun
  # Strings
  ps.ids = ParamHelpers::getParamIds(ps.mbo, repeated = TRUE, with.nr = TRUE)
  y.name.mbo = ctrl.mbo$y.name
  infill.mbo = ctrl.mbo$infill.crit$id
  if(infill.mbo == "cb") cb.lambda = ctrl.mbo$infill.crit$params$cb.lambda
  # Numeric
  # stored are the stored models in the process, not every iteration has a stored model
  stored = sort(as.integer(names(res.mbo$models)))
  iters.mbo = opdf[nrow(opdf), "dob"] #max(ParamHelpers::getOptPathDOB(op.mbo)),op.mbo = res.mbo$opt.path
  init.size = nrow(res.mbo$final.opt.state$opt.problem$design)
  maximize.mult = if (ctrl.mbo$minimize) 1 else -1 # needed to construct phi(cb) in mergeShapleyRes
  # Tables
  if(infill.mbo == "cb"){
    props = opdf[which(opdf$dob != 0), c(ps.ids, infill.mbo, "se", "mean"), drop = FALSE]
  } else {
    props = opdf[which(opdf$dob != 0), c(ps.ids, infill.mbo), drop = FALSE]
  }
  rownames(props) = NULL
  
  # Input checking
  # stops if multi point proposal or multiobjective target fun or max. problems
  if (ctrl.mbo$propose.points > 1 || ctrl.mbo$n.objectives > 1) 
    stop("ShapleyMBO not implemented for Multipoint-Proposal or Multiobjective function")
  # stops for maximization problems - not yet robust for such cases
  if (ctrl.mbo$minimize == FALSE) 
    stop("ShapleyMBO not implemented and tested yet for Max problmes")
  # stops if infill is not one of the seven built-in infill in mlrMBO
  if (!(infill.mbo %in% c("mean", "se", "ei", "cb", "eqi", "aei", "adacb"))) 
    stop("ShapleyMBO only implemented for the seven built-in single obj. infill crits")
  if (infill.mbo != "cb")
    warning("ShapleyMBO implemented bit not yet tested in detail for infill criteria other than cb")
  
  # iter.interest
  checkmate::assertNumeric(iter.interest, lower = 1, upper = iters.mbo, any.missing = FALSE,
                           all.missing = FALSE, min.len = 1, max.len = iters.mbo, unique = TRUE,
                           null.ok = TRUE)
  #iter.interest need to be a subset of stored
  checkmate::assertSubset(iter.interest, stored, empty.ok = TRUE)
  
  # stored surrogate models 
  if (length(stored) > 1 && max(stored) == iters.mbo + 1) {
    # we need to remove the final SM, since there are no further BO steps
    stored = stored[-length(stored)]
  }
  if(is.null(iter.interest)) stored = stored # if NULL analyze all possible iterations
  if(!(is.null(iter.interest))) stored = sort(stored[iter.interest]) # otherwise only selected ones
  
  # contribution - decomposition can be applied only with the CB
  checkmate::assertLogical(contribution, len = 1)
  if(contribution == TRUE && infill.mbo != "cb") 
    stop("Shapley Value decomposition only available for the Confidence Bound infill")
  # seed
  checkmate::assertNumber(seed)
  # no.cores
  checkmate::assertNumber(no.cores)
  
  ## Begin of the Computation ##
  if(contribution == FALSE ) {
    res = parallel::mclapply(
      stored,
      function(x) {
        set.seed(seed)
        ShapleyAf(
          res.mbo.p = res.mbo, # the results of the mbo run are required because for each iter we need to create the PredictorAf object
          iter.p = x,
          x.interest.s = props[x, c(ps.ids, infill.mbo), drop = FALSE],
          sample.size.s = sample.size
        )
      },
      mc.cores = no.cores
    )
    names(res) = stored
    res = dplyr::bind_rows(res, .id = "iter") %>% as.data.frame()
  }
  
  if(contribution == TRUE && infill.mbo == "cb") {
    # MEAN contributions
    res.mean = parallel::mclapply(
      stored,
      function(x) {
        set.seed(seed)
        ShapleyMean(
          model.p = sm[[as.character(x)]],
          x.interest.s = props[x, c(ps.ids, "mean"), drop = FALSE],
          sample.size.s = sample.size, 
          ps.mbo = ps.mbo
        )
      },
      mc.cores = no.cores
    )
    # each list element is named with the corresponding iteration
    names(res.mean) = stored
    
    # SE contributions
    # repeat the same procedure
    res.se = parallel::mclapply(
      stored,
      function(x) {
        set.seed(seed)
        ShapleySe(
          model.p = sm[[as.character(x)]],
          x.interest.s = props[x, c(ps.ids, "se"), drop = FALSE],
          sample.size.s = sample.size,
          ps.mbo = ps.mbo
        )
      },
      mc.cores = no.cores
    )
    # each list element is named with the corresponding iteration
    names(res.se) = stored
    
    # merge the Mean and Se results together and compute the SV of the CB
    res = mergeShapleyRes(res.mean, res.se, 
                         lambda = cb.lambda,
                         max.mult = maximize.mult, 
                         sample.size.s = sample.size
          )
    res = as.data.frame(res)
  }
  
  return(res)
}
