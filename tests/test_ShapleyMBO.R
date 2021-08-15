######################### ShapleyMBO ########################################
######################### TEST ##############################################
library(testthat)
library(mlrMBO)
library(iml)

# objective function
set.seed(1)
obj.fun = makeHyperEllipsoidFunction(dimensions = 4L)
ps = getParamSet(obj.fun)
des = generateDesign(n = 16, par.set = ps, fun = lhs::maximinLHS)
ctrl = makeMBOControl(y.name = "y", store.model.at = 1:4)
ctrl = setMBOControlInfill(ctrl, opt = "focussearch", crit = makeMBOInfillCritCB(cb.lambda = 10))
ctrl = setMBOControlTermination(ctrl, iters = 3)
res = mbo(obj.fun, design = des, control = ctrl, show.info = FALSE)

test_that("ShapleyMBO is reproducible", {
  # no contribution
  S1 = ShapleyMBO(res.mbo = res, seed = 123)
  S2 = ShapleyMBO(res.mbo = res, seed = 123)
  expect_identical(S1, S2)

  # with contribution
  S3 = ShapleyMBO(res.mbo = res, contribution = TRUE, seed = 123)
  S4 = ShapleyMBO(res.mbo = res, contribution = TRUE, seed = 123)
  expect_identical(S3, S4)
  
  # with contr, default seed
  S5 = ShapleyMBO(res.mbo = res, contribution = TRUE)
  S6 = ShapleyMBO(res.mbo = res, contribution = TRUE)
  expect_identical(S5, S6)
  
  # w/o contr, default seed
  S5 = ShapleyMBO(res.mbo = res, contribution = FALSE)
  S6 = ShapleyMBO(res.mbo = res, contribution = FALSE)
  expect_identical(S5, S6)
  
})

test_that("stored and assertions on iter.interest in ShapleyMBO works correctly", {
  S1 = ShapleyMBO(res.mbo = res, seed = 1, iter.interest = 2)
  expect_equal(
    unique(S1$iter), "2"
  )
  # iters are sorted
  S2 = ShapleyMBO(res.mbo = res, seed = 1, iter.interest = c(3,1,2))
  expect_equal(
    unique(S2$iter), c("1", "2", "3")
  )
  # last surrogate model is removed from stored
  S3 = ShapleyMBO(res.mbo = res, seed = 1, iter.interest = NULL)
  expect_equal(
    unique(S3$iter), c("1", "2", "3")
  )
  
  expect_error(ShapleyMBO(res.mbo = res, seed = 1, iter.interest = c(1,5)))
  expect_error(ShapleyMBO(res.mbo = res, seed = 1, iter.interest = c(3,3,3)))

  expect_error(
    ShapleyMBO(res.mbo = res, seed = 1, iter.interest = 1:5)
  )
})

test_that("Linearity Axiom works - cb: phi, pred.interest, pred.average", {
  # test e.g. for iter 2
  S = ShapleyMBO(res.mbo = res, iter.interest = 2, contribution = FALSE, seed = 123)
  opdf = as.data.frame(res$opt.path)
  x_mean = opdf[18, c("x1", "x2", "x3", "x4", "mean"), drop = FALSE]
  sm = res$models$`2`
  # compute the Shapley Value for the mean
  set.seed(123)
  samples = ParamHelpers::generateDesign(n = 4000, par.set = res$opt.path$par.set, fun = lhs::randomLHS)
  data_mean = cbind(samples, mean = predict(sm, newdata = samples)$data$response)
  X_mean = rbind(
    x_mean,
    data_mean,
    make.row.names = FALSE
  )
  P_mean = iml::Predictor$new(model = sm, data = X_mean, y = "mean")
  S_mean = iml::Shapley$new(predictor = P_mean, 
                            x.interest = X_mean[1, ], 
                            sample.size = 100)
  res_mean = S_mean$results %>% dplyr::mutate(
    pred.interest = S_mean$y.hat.interest,
    pred.average = S_mean$y.hat.average
  )
  
  # compute the Shapley Value for the se
  x_se = opdf[18, c("x1", "x2", "x3", "x4", "se"), drop = FALSE]
  set.seed(123)
  samples = ParamHelpers::generateDesign(n = 4000, par.set = res$opt.path$par.set, fun = lhs::randomLHS)
  data_se = cbind(samples, se = predict(sm, newdata = samples)$data$se)
  X_se = rbind(
    x_se,
    data_se,
    make.row.names = FALSE
  )
  predfun = function(model, newdata) {
    .prediction = predict(sm, newdata = newdata)$data$se
  }
  P_se = iml::Predictor$new(model = NULL, predict.function = predfun, data = X_se, y = "se")
  S_se = iml::Shapley$new(predictor = P_se, 
                          x.interest = X_se[1, ], 
                          sample.size = 100)
  res_se = S_se$results %>% dplyr::mutate(
    pred.interest = S_se$y.hat.interest,
    pred.average = S_se$y.hat.average
  )
  
  lambda = res$control$infill.crit$params$cb.lambda
  mm = ifelse(res$control$minimize, 1, -1)
  
  expect_equal(S[,"pred.interest"], mm * res_mean[, "pred.interest"] - lambda * res_se[,"pred.interest"])
  expect_equal(S[,"pred.average"], mm * res_mean[, "pred.average"] - lambda * res_se[,"pred.average"])
  expect_equal(S[,"phi"], mm * res_mean[, "phi"] - lambda * res_se[,"phi"])
  
})

test_that("Shapley objects have the same design for mean and se WITHIN iters", {
  S = ShapleyMBO_env(res.mbo = res, contribution = TRUE, seed = 123)
  # iter1
  expect_identical(
    S$mean$`1`$.__enclos_env__$private$dataDesign,
    S$se$`1`$.__enclos_env__$private$dataDesign
  )
  expect_identical(
    S$mean$`1`$.__enclos_env__$private$dataSample,
    S$se$`1`$.__enclos_env__$private$dataSample
  )
  #iter2
  expect_identical(
    S$mean$`2`$.__enclos_env__$private$dataDesign,
    S$se$`2`$.__enclos_env__$private$dataDesign
  )
  expect_identical(
    S$mean$`2`$.__enclos_env__$private$dataSample,
    S$se$`2`$.__enclos_env__$private$dataSample
  )
  #iter3
  expect_identical(
    S$mean$`3`$.__enclos_env__$private$dataDesign,
    S$se$`3`$.__enclos_env__$private$dataDesign
  )
  expect_identical(
    S$mean$`3`$.__enclos_env__$private$dataSample,
    S$se$`3`$.__enclos_env__$private$dataSample
  )
})

test_that("Shapley objects have the setting for mean and se BETWEEN iters", {
  S = ShapleyMBO_env(res.mbo = res, contribution = TRUE, seed = 123)
  # here we can only test dataSample since as we explain different instances in each iteration
  # dataDesign will differ nevertheless between iterations
  # mean
  #dataSample
  expect_identical(
      S$mean$`1`$.__enclos_env__$private$dataSample,
      S$mean$`2`$.__enclos_env__$private$dataSample
  )
  expect_identical(
      S$mean$`1`$.__enclos_env__$private$dataSample,
      S$mean$`3`$.__enclos_env__$private$dataSample
  )
  expect_identical(
      S$mean$`2`$.__enclos_env__$private$dataSample,
      S$mean$`3`$.__enclos_env__$private$dataSample
  )

  # se
  #dataSample
  expect_identical(
      S$se$`1`$.__enclos_env__$private$dataSample,
      S$se$`2`$.__enclos_env__$private$dataSample
  )
  expect_identical(
      S$se$`1`$.__enclos_env__$private$dataSample,
      S$se$`3`$.__enclos_env__$private$dataSample
  )
  expect_identical(
      S$se$`2`$.__enclos_env__$private$dataSample,
      S$se$`3`$.__enclos_env__$private$dataSample
  )
  contribution = FALSE
  S2 = ShapleyMBO_env(res.mbo = res, contribution = FALSE, seed = 123)
  expect_identical(
    S2$`1`$.__enclos_env__$private$dataSample,
    S2$`2`$.__enclos_env__$private$dataSample
  )
  expect_identical(
    S2$`1`$.__enclos_env__$private$dataSample,
    S2$`3`$.__enclos_env__$private$dataSample
  )
  expect_identical(
    S2$`2`$.__enclos_env__$private$dataSample,
    S2$`3`$.__enclos_env__$private$dataSample
  )
})

test_that("utils_ShapleyMBO correct - ShapleyMBO gives same results with contribution TRUE and FALSE", {
  
  S1 = ShapleyMBO(res.mbo = res, seed = 123, contribution = FALSE)
  S2 = ShapleyMBO(res.mbo = res, seed = 123, contribution = TRUE)
  df = S2[, c(1,2,16,17,3,14,15)] # select the interesting columns
  colnames(df) = colnames(S1)
  expect_equal(S1, df)

  #w/o seed
  S3 = ShapleyMBO(res.mbo = res, contribution = FALSE)
  S4 = ShapleyMBO(res.mbo = res, contribution = TRUE)
  df = S4[, c(1,2,16,17,3,14,15)]
  colnames(df) = colnames(S3)
  expect_equal(S3, df)
})

test_that("ShapleyMBO explains the correct instances", {
  # contribution = TRUE
  S = ShapleyMBO(res.mbo = res, seed = 123, contribution = TRUE)
  S = S %>% tidyr::separate(feature.value, into = c("f", "feature.value"), sep ="=", convert = TRUE)
  explained1 = c(S[1:4, "feature.value"], 
                 mean = S[1, "pred.interest_mean"], 
                 se = S[1, "pred.interest_se"], 
                 cb = S[1, "pred.interest_cb"]
  )
  names(explained1)[1:4] = c("x1", "x2", "x3", "x4")
  explained2 = c(S[5:8, "feature.value"], 
                 mean = S[5, "pred.interest_mean"], 
                 se = S[5, "pred.interest_se"], 
                 cb = S[5, "pred.interest_cb"]
  )
  names(explained2)[1:4] = c("x1", "x2", "x3", "x4")
  explained3 = c(S[9:12, "feature.value"], 
                 mean = S[9, "pred.interest_mean"], 
                 se = S[9, "pred.interest_se"], 
                 cb = S[9, "pred.interest_cb"]
  )
  names(explained3)[1:4] = c("x1", "x2", "x3", "x4")
  
  opdf = as.data.frame(res$opt.path) %>% dplyr::select(x1, x2, x3, x4, mean, se, cb)
  x.interest1 = opdf[17, , drop = TRUE] %>% unlist()
  x.interest2 = opdf[18, , drop = TRUE] %>% unlist()
  x.interest3 = opdf[19, , drop = TRUE] %>% unlist()
  
  expect_equal(explained1, x.interest1)
  expect_equal(explained2, x.interest2)
  expect_equal(explained3, x.interest3)
})

test_that("ShapleyMBO and ShapleyMBO_mclapply give the same results", {
  S1 = ShapleyMBO(res.mbo = res, contribution = TRUE)
  S2 = ShapleyMBO_mclapply(res.mbo = res, contribution = TRUE, no.cores = 3)
  expect_identical(S1, S2)

  S3 = ShapleyMBO(res.mbo = res, contribution = TRUE, seed = 123)
  S4 = ShapleyMBO_mclapply(res.mbo = res, contribution = TRUE, no.cores = 3, seed = 123)
  expect_identical(S3, S4)

  S5 = ShapleyMBO(res.mbo = res, contribution = FALSE)
  S6 = ShapleyMBO_mclapply(res.mbo = res, contribution = FALSE, no.cores = 3)
  expect_identical(S5, S6)
})
#- mm.mbo works correctly



