###################################################### 
############       ShapleySe     #####################
######################################################
#- computes the Shapley Values of the parameter of a MBO proposal in a given iteration, 
#  using the predicted model uncertainty ($se) as payout
#- this is equivalent to a standard iml procedure with custom prediction function
ShapleySe = function(model.p,
                     x.interest.s, # explicand, instance to explain 
                     sample.size.s = 100,
                     ps.mbo
) {
  # 1. sample instances from the input space
  ps.ids = ParamHelpers::getParamIds(ps.mbo, repeated = TRUE, with.nr = TRUE)
  N = 1000 * length(ps.ids)
  data.p = ParamHelpers::generateDesign(n = N, par.set = ps.mbo, fun = lhs::randomLHS)
  
  # 2. create the data to pass to the Predictor object
  data = cbind(data.p, se = predict(model.p, newdata = data.p)$data$se)
  X = rbind(
    x.interest.s, # explicand is the first row
    data,
    make.row.names = FALSE
  )
  # 3. create Predictor Object
  # for the se we create a custom predict function
  predfun = function(model, newdata) {
    .prediction = predict(model.p, newdata = newdata)$data$se
  }
  P = iml::Predictor$new(model = NULL, predict.function = predfun, data = X, y = "se")
  
  # 4. create the Shapley Object
  S = iml::Shapley$new(
    predictor = P, x.interest = X[1, ], sample.size = sample.size.s
    )
  
  # 5. get the results
  res = getShapleyRes(S, "detailed")
  
  return(res)
}
