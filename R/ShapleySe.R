###################################################### 
############       ShapleySe     #####################
######################################################
#- computes the Shapley Values of the feature for MBO proposal in a given iteration, 
#  using the model uncertainty as payout
#- this is like a standard iml procedure with custom prediction function
ShapleySe = function(model.p,
                     x.interest.s, # explicand, instance to explain 
                     sample.size.s = 100,
                     ps.mbo
) {
  # 1. sample instances from the input space, passed to the Predictor later on
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
  # if(is.null(seed)) set.seed(1)
  # else set.seed(seed) #seed is needed is order to use the same coalitions and feature values between mean and se
  
  S = iml::Shapley$new(
    predictor = P, x.interest = X[1, ], sample.size = sample.size.s
    )
  
  # 5. get the results
  res = getShapleyRes(S, "detailed")
  
  return(res)
}
