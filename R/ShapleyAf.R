###################################################### 
##########       ShapleyAf       #####################
######################################################
#- computes the Shapley Values of the parameters of the MBO proposal in a given iteration, 
#  using the af as payout
ShapleyAf = function(res.mbo.p, 
                     iter.p, 
                     x.interest.s, # explicand, here we only need the infill column
                     sample.size.s = 100
) {
  
  # extract information
  opdf = as.data.frame(res.mbo.p$opt.path)
  # Par.Set and Ctrl Object
  ps = res.mbo.p$opt.path$par.set
  ctrl = res.mbo.p$control
  # Pred Models
  sm = res.mbo.p$models[[as.character(iter.p)]]
  af = ctrl$infill.crit$fun
  # Strings
  ps.ids = ParamHelpers::getParamIds(ps, repeated = TRUE, with.nr = TRUE)
  y.name.mbo = ctrl$y.name
  infill = ctrl$infill.crit$id
  # Other
  init.size = nrow(res.mbo.p$final.opt.state$opt.problem$design)
  
  # 1. Compute the design
  des = opdf[1:(init.size + iter.p - 1), c(ps.ids, y.name.mbo), drop = FALSE]
  
  # 2. sample instances from the input space
  N = 1000 * length(ps.ids)
  data.p = ParamHelpers::generateDesign(n = N, par.set = ps, fun = lhs::randomLHS)
  
  # 3. compute AF value for the sampled instances
  data.p$af = af(
    points = data.p,
    models = list(sm),
    control = ctrl,
    par.set = ps,
    designs = des,
    iter = iter.p, 
    progress = getProgressAdaCB(res.mbo = res.mbo.p, iter = iter.p),
    attributes = FALSE
  )
  colnames(data.p)[colnames(data.p) == "af"] = infill
  
  # 4. bind the proposal (in the first row) and the samples
  X = rbind(x.interest.s,
            data.p,
            make.row.names = FALSE
  )
  
  # 5. create a Predictor object
  P = PredictorAf$new(
    model = sm, data = X, y = infill, #batch.size = batch.size.p,
    res.mbo = res.mbo.p, design = des, iter = iter.p
  )
  # 6. compute the Shapley Values
  S = iml::Shapley$new(
    predictor = P, x.interest = X[1, ], sample.size = sample.size.s
  )
  
  # 7. get the results
  res = getShapleyRes(S, "simple")
  
  return(res)
}
