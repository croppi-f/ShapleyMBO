###################################################### 
##########       ShapleyMean     #####################
######################################################
#- computes the Shapley Values of the paarmeter of MBO proposal in a given iteration, 
#  using the model prediction  ($response) as payout
#- this is equivalent to the standard iml procedure
ShapleyMean = function(model.p,
                       x.interest.s, # explicand, instance to explain 
                       sample.size.s = 100,
                       ps.mbo
) {
  # 1. sample instances from the input space
  ps.ids = ParamHelpers::getParamIds(ps.mbo, repeated = TRUE, with.nr = TRUE)
  N = 1000 * length(ps.ids)
  data.p = ParamHelpers::generateDesign(n = N, par.set = ps.mbo, fun = lhs::randomLHS)
  
  # 2. create the data to pass to the Predictor object
  data = cbind(data.p, mean = predict(model.p, newdata = data.p)$data$response)
  X = rbind(
    x.interest.s, # explicand is the first row
    data,
    make.row.names = FALSE
  )
  # 3. create Predictor Object
  P = iml::Predictor$new(model = model.p, data = X, y = "mean")
  
  # 4. create the Shapley Object
  S = iml::Shapley$new(
    predictor = P, x.interest = X[1, ], sample.size = sample.size.s
    )
  
  # 5. get the results
  res = getShapleyRes(S, "detailed")
  
  return(res)
}
