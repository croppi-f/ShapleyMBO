################# Data and Configurations for the tests ##############################
# source the data
source("tests/test_data/md.R")
library(mlrMBO)
model = train(learner = makeLearner("regr.randomForest"), 
              task = makeRegrTask(data = data, target = "target")
)

# Bayesian Optimization
# getPredictionResponse are the evaliation of the (new) proposed points
fun = function(x) {
  df = as.data.frame(x)
  return(getPredictionResponse(predict(model, newdata = df)))
}

ps = makeParamSet(
  makeNumericParam("f", lower = 0, upper = 0.25),
  makeNumericParam("k", lower = 0, upper = 0.1),
  makeNumericParam("du", lower = 0, upper = 2e-5),
  makeNumericParam("dv", lower = 0, upper = 2e-5),
  makeIntegerParam("x", lower = 0, upper = 200),
  makeIntegerParam("y", lower = 0, upper = 200)
)

objfun = makeSingleObjectiveFunction(
  name = "Synthesis",
  fn = fun,
  par.set = ps,
  has.simple.signature = FALSE,
  minimize = TRUE
)

# sample 10 points for the initial surrogate model
# initial.data = data[sample(1:nrow(data), 10), ]
set.seed(1)
initial.data = generateDesign(10, attributes(objfun)$par.set)

#source modified function from mlrMBO
source("R/utils_xplxpl/proposePointsByInfillOptimization.R")
source("R/utils_xplxpl/makeMBOResult.OptState.R")
source("R/utils_xplxpl/getSupportedInfillOptFunctions.R")
source("R/utils_xplxpl/proposePointsByInfillOptimization.R")
source("R/utils_xplxpl/getInfillOptFunction.R")
source("R/utils_xplxpl/checkStuff.R")

# source new infill optimization functions "...Savepts"
source("R/utils_xplxpl/infillOptFocusSavepts.R")
source("R/utils_xplxpl/infillOptEASavepts.R")
source("R/utils_xplxpl/infillOptCMAESSavepts.R")

######################### TEST #################################################
library(testthat)

test_that("Equivalence for Confidence Bound", {
  ctrl = makeMBOControl(y.name = "target", store.model.at = 1:3)
  ctrl = setMBOControlTermination(ctrl, iters = 2)
  ctrl = setMBOControlInfill(ctrl,
                                opt = "focussearchSavepts",
                                opt.focussearch.maxit = 20,
                                opt.focussearch.points = 50,
                                crit = makeMBOInfillCritCB()
  )
  res = mbo(objfun,
                   design = initial.data,
                   control = ctrl,
                   show.info = FALSE
  )
  opdf = as.data.frame(res$opt.path)
  pars = names(res$x)
  y.name = res$control$y.name
  infill = res$control$infill.crit$id
  pp.true = opdf[which(opdf$dob > 0), c(pars, y.name)]
  
  #create the PredictorAfObject in iter 1
  p1 = PredictorAf$new(model = res$models[[1]],
                         data = res$seen.points[[1]], 
                         res.mbo = res,
                         iter = 1,
                         design = opdf[1:10, c(pars, y.name)] 
  )
  #create the PredictorAfObject in iter 2
  p2 = PredictorAf$new(model = res$models[[2]],
                       data = res$seen.points[[2]], 
                       res.mbo = res,
                       iter = 2,
                       design = opdf[1:11, c(pars, y.name)] 
  )
  # iter 1
  expect_equivalent(
    p1$predict(pp.true[1, ]),
    opdf[11, infill]
  )
  # iter 2
  expect_equivalent(
    p2$predict(pp.true[2, ]),
    opdf[12, infill]
  )
})

test_that("Equivalence for Expected Improvement", {
  ctrl = makeMBOControl(y.name = "target", store.model.at = 1:3)
  ctrl = setMBOControlTermination(ctrl, iters = 2)
  ctrl = setMBOControlInfill(ctrl,
                             opt = "focussearchSavepts",
                             opt.focussearch.maxit = 20,
                             opt.focussearch.points = 50,
                             crit = makeMBOInfillCritEI()
  )
  res = mbo(objfun,
            design = initial.data,
            control = ctrl,
            show.info = FALSE
  )
  opdf = as.data.frame(res$opt.path)
  pars = names(res$x)
  y.name = res$control$y.name
  infill = res$control$infill.crit$id
  pp.true = opdf[which(opdf$dob > 0), c(pars, y.name)]
  
  #create the PredictorAfObject in iter 1
  p1 = PredictorAf$new(model = res$models[[1]],
                       data = res$seen.points[[1]], 
                       res.mbo = res,
                       iter = 1,
                       design = opdf[1:10, c(pars, y.name)] 
  )
  #create the PredictorAfObject in iter 2
  p2 = PredictorAf$new(model = res$models[[2]],
                       data = res$seen.points[[2]], 
                       res.mbo = res,
                       iter = 2,
                       design = opdf[1:11, c(pars, y.name)] 
  )
  # iter 1
  expect_equivalent(
    p1$predict(pp.true[1, ]),
    opdf[11, infill]
  )
  # iter 2
  expect_equivalent(
    p2$predict(pp.true[2, ]),
    opdf[12, infill]
  )
  
})

test_that("Equivalence for Mean Response", {
  ctrl = makeMBOControl(y.name = "target", store.model.at = 1:3)
  ctrl = setMBOControlTermination(ctrl, iters = 2)
  ctrl = setMBOControlInfill(ctrl,
                             opt = "focussearchSavepts",
                             opt.focussearch.maxit = 20,
                             opt.focussearch.points = 50,
                             crit = makeMBOInfillCritMeanResponse()
  )
  res = mbo(objfun,
            design = initial.data,
            control = ctrl,
            show.info = FALSE
  )
  opdf = as.data.frame(res$opt.path)
  pars = names(res$x)
  y.name = res$control$y.name
  infill = res$control$infill.crit$id
  pp.true = opdf[which(opdf$dob > 0), c(pars, y.name)]
  
  #create the PredictorAfObject in iter 1
  p1 = PredictorAf$new(model = res$models[[1]],
                       data = res$seen.points[[1]], 
                       res.mbo = res,
                       iter = 1,
                       design = opdf[1:10, c(pars, y.name)] 
  )
  #create the PredictorAfObject in iter 2
  p2 = PredictorAf$new(model = res$models[[2]],
                       data = res$seen.points[[2]], 
                       res.mbo = res,
                       iter = 2,
                       design = opdf[1:11, c(pars, y.name)] 
  )
  # iter 1
  expect_equivalent(
    p1$predict(pp.true[1, ]),
    opdf[11, infill]
  )
  # iter 2
  expect_equivalent(
    p2$predict(pp.true[2, ]),
    opdf[12, infill]
  )
})

test_that("Equivalence for Standard Error", {
  ctrl = makeMBOControl(y.name = "target", store.model.at = 1:3)
  ctrl = setMBOControlTermination(ctrl, iters = 2)
  ctrl = setMBOControlInfill(ctrl,
                             opt = "focussearchSavepts",
                             opt.focussearch.maxit = 20,
                             opt.focussearch.points = 50,
                             crit = makeMBOInfillCritStandardError()
  )
  res = mbo(objfun,
            design = initial.data,
            control = ctrl,
            show.info = FALSE
  )
  opdf = as.data.frame(res$opt.path)
  pars = names(res$x)
  y.name = res$control$y.name
  infill = res$control$infill.crit$id
  pp.true = opdf[which(opdf$dob > 0), c(pars, y.name)]
  
  #create the PredictorAfObject in iter 1
  p1 = PredictorAf$new(model = res$models[[1]],
                       data = res$seen.points[[1]], 
                       res.mbo = res,
                       iter = 1,
                       design = opdf[1:10, c(pars, y.name)] 
  )
  #create the PredictorAfObject in iter 2
  p2 = PredictorAf$new(model = res$models[[2]],
                       data = res$seen.points[[2]], 
                       res.mbo = res,
                       iter = 2,
                       design = opdf[1:11, c(pars, y.name)] 
  )
  # iter 1
  expect_equivalent(
    p1$predict(pp.true[1, ]),
    opdf[11, infill]
  )
  # iter 2
  expect_equivalent(
    p2$predict(pp.true[2, ]),
    opdf[12, infill]
  )
})

test_that("Equivalence for Augmented EI", {
  ctrl = makeMBOControl(y.name = "target", store.model.at = 1:3)
  ctrl = setMBOControlTermination(ctrl, iters = 2)
  ctrl = setMBOControlInfill(ctrl,
                             opt = "focussearchSavepts",
                             opt.focussearch.maxit = 20,
                             opt.focussearch.points = 50,
                             crit = makeMBOInfillCritAEI()
  )
  res = mbo(objfun,
            design = initial.data,
            control = ctrl,
            show.info = FALSE
  )
  opdf = as.data.frame(res$opt.path)
  pars = names(res$x)
  y.name = res$control$y.name
  infill = res$control$infill.crit$id
  pp.true = opdf[which(opdf$dob > 0), c(pars, y.name)]
  
  #create the PredictorAfObject in iter 1
  p1 = PredictorAf$new(model = res$models[[1]],
                       data = res$seen.points[[1]], 
                       res.mbo = res,
                       iter = 1,
                       design = opdf[1:10, c(pars, y.name)] 
  )
  #create the PredictorAfObject in iter 2
  p2 = PredictorAf$new(model = res$models[[2]],
                       data = res$seen.points[[2]], 
                       res.mbo = res,
                       iter = 2,
                       design = opdf[1:11, c(pars, y.name)] 
  )
  # iter 1
  expect_equivalent(
    p1$predict(pp.true[1, ]),
    opdf[11, infill]
  )
  # iter 2
  expect_equivalent(
    p2$predict(pp.true[2, ]),
    opdf[12, infill]
  )
})

test_that("Equivalence for Expected Quantile Improvement", {
  ctrl = makeMBOControl(y.name = "target", store.model.at = 1:3)
  ctrl = setMBOControlTermination(ctrl, iters = 2)
  ctrl = setMBOControlInfill(ctrl,
                             opt = "focussearchSavepts",
                             opt.focussearch.maxit = 20,
                             opt.focussearch.points = 50,
                             crit = makeMBOInfillCritEQI()
  )
  res = mbo(objfun,
            design = initial.data,
            control = ctrl,
            show.info = FALSE
  )
  opdf = as.data.frame(res$opt.path)
  pars = names(res$x)
  y.name = res$control$y.name
  infill = res$control$infill.crit$id
  pp.true = opdf[which(opdf$dob > 0), c(pars, y.name)]
  
  #create the PredictorAfObject in iter 1
  p1 = PredictorAf$new(model = res$models[[1]],
                       data = res$seen.points[[1]], 
                       res.mbo = res,
                       iter = 1,
                       design = opdf[1:10, c(pars, y.name)] 
  )
  #create the PredictorAfObject in iter 2
  p2 = PredictorAf$new(model = res$models[[2]],
                       data = res$seen.points[[2]], 
                       res.mbo = res,
                       iter = 2,
                       design = opdf[1:11, c(pars, y.name)] 
  )
  # iter 1
  expect_equivalent(
    p1$predict(pp.true[1, ]),
    opdf[11, infill]
  )
  # iter 2
  expect_equivalent(
    p2$predict(pp.true[2, ]),
    opdf[12, infill]
  )
})

test_that("Equivalence for Adaptive CB & term.cond. = iters", {
  ctrl = makeMBOControl(y.name = "target", store.model.at = 1:3)
  ctrl = setMBOControlTermination(ctrl, iters = 2)
  ctrl = setMBOControlInfill(ctrl,
                             opt = "focussearchSavepts",
                             opt.focussearch.maxit = 20,
                             opt.focussearch.points = 50,
                             crit = makeMBOInfillCritEI()
  )
  res = mbo(objfun,
            design = initial.data,
            control = ctrl,
            show.info = FALSE
  )
  opdf = as.data.frame(res$opt.path)
  pars = names(res$x)
  y.name = res$control$y.name
  infill = res$control$infill.crit$id
  pp.true = opdf[which(opdf$dob > 0), c(pars, y.name)]
  
  #create the PredictorAfObject in iter 1
  p1 = PredictorAf$new(model = res$models[[1]],
                       data = res$seen.points[[1]], 
                       res.mbo = res,
                       iter = 1,
                       design = opdf[1:10, c(pars, y.name)] 
  )
  #create the PredictorAfObject in iter 2
  p2 = PredictorAf$new(model = res$models[[2]],
                       data = res$seen.points[[2]], 
                       res.mbo = res,
                       iter = 2,
                       design = opdf[1:11, c(pars, y.name)] 
  )
  # iter 1
  expect_equivalent(
    p1$predict(pp.true[1, ]),
    opdf[11, infill]
  )
  # iter 2
  expect_equivalent(
    p2$predict(pp.true[2, ]),
    opdf[12, infill]
  )
})
## FIXME: 
test_that("Equivalence for AdaCB over multiple iterations and term.cond. time budget", {
  # create the PredictorAfObject in iter 1
  p.adacb.tb = PredictorAf$new(model = res.mbo.adacb.tb$models[[1]], 
                                 data = res.mbo.adacb.tb$seen.points[[1]], 
                                 res.mbo = res.mbo.adacb.tb,
                                 iter = 1, 
                                 design = opdf.adacb.tb[1:20,c(names(res.mbo.adacb.tb$x), res.mbo.adacb.tb$control$y.name)]
  )
  # iter 3
  p3.adacb.tb = PredictorAf$new(model = res.mbo.adacb.tb$models[[3]], 
                               data = res.mbo.adacb.tb$seen.points[[3]], 
                               res.mbo = res.mbo.adacb.tb,
                               iter = 3, 
                               design = opdf.adacb.tb[1:22,c(names(res.mbo.adacb.tb$x), res.mbo.adacb.tb$control$y.name)]
  )
  
  # iter 7
  p7.adacb.tb = PredictorAf$new(model = res.mbo.adacb.tb$models[[7]], 
                                data = res.mbo.adacb.tb$seen.points[[7]], 
                                res.mbo = res.mbo.adacb.tb,
                                iter = 7, 
                                design = opdf.adacb.tb[1:26,c(names(res.mbo.adacb.tb$x), res.mbo.adacb.tb$control$y.name)]
  )
  
  expect_equivalent(
    p.adacb.tb$predict(pp.true.adacb.tb[1,1:6]),
    opdf.adacb.tb[21, res.mbo.adacb.tb$control$infill.crit$id]
  )
  expect_equivalent(
    p3.adacb.tb$predict(pp.true.adacb.tb[3,1:6]),
    opdf.adacb.tb[23, res.mbo.adacb.tb$control$infill.crit$id]
  )
  expect_equivalent(
    p7.adacb.tb$predict(pp.true.adacb.tb[7,1:6]),
    opdf.adacb.tb[27, res.mbo.adacb.tb$control$infill.crit$id]
  )
  
})

test_that("Data corresponds to seen points and designs to training set within a specific iteration", {
  ctrl = makeMBOControl(y.name = "target", store.model.at = 1:3)
  ctrl = setMBOControlTermination(ctrl, iters = 2)
  ctrl = setMBOControlInfill(ctrl,
                             opt = "focussearchSavepts",
                             opt.focussearch.maxit = 20,
                             opt.focussearch.points = 50,
                             crit = makeMBOInfillCritEI()
  )
  res = mbo(objfun,
            design = initial.data,
            control = ctrl,
            show.info = FALSE
  )
  opdf = as.data.frame(res$opt.path)
  pars = names(res$x)
  y.name = res$control$y.name
  infill = res$control$infill.crit$id
  pp.true = opdf[which(opdf$dob > 0), c(pars, y.name)]
  # Iter 1
  #create the PredictorAfObject in iter 1
  p1 = PredictorAf$new(model = res$models[[1]],
                       data = res$seen.points[[1]], 
                       res.mbo = res,
                       iter = 1,
                       design = opdf[1:10, c(pars, y.name)] 
  )
  
  expect_equal(
    p1$data$X, data.table::data.table(res$seen.points[[1]])
  )
  expect_equal(
    p1$design[[1]], opdf[1:10,c(pars, y.name)]
  )
  
  # Iter 2
  p2 = PredictorAf$new(model = res$models[[2]],
                       data = res$seen.points[[2]], 
                       res.mbo = res,
                       iter = 2,
                       design = opdf[1:11, c(pars, y.name)] 
  )
  
  expect_equal(
    p2$data$X, data.table::data.table(res$seen.points[[2]])
  )
  expect_equal(
    p2$design[[1]], opdf[1:11,c(pars, y.name)]
  )
})

test_that("PredictorAf works also if y not specified in data or as an argument", {
  
  ctrl = makeMBOControl(y.name = "target", store.model.at = 1:3)
  ctrl = setMBOControlTermination(ctrl, iters = 2)
  ctrl = setMBOControlInfill(ctrl,
                             opt = "focussearchSavepts",
                             opt.focussearch.maxit = 20,
                             opt.focussearch.points = 50,
                             crit = makeMBOInfillCritEI()
  )
  res = mbo(objfun,
            design = initial.data,
            control = ctrl,
            show.info = FALSE
  )
  
  opdf = as.data.frame(res$opt.path)
  pars = names(res$x)
  y.name = res$control$y.name
  infill = res$control$infill.crit$id
  pp.true = opdf[which(opdf$dob > 0), c(pars, y.name)]
  acq.fun = res$control$infill.crit$fun
  
  # iter 1 
  # without y
  p1 = PredictorAf$new(model = res$models[[1]],
                       data = res$seen.points[[1]], 
                       res.mbo = res,
                       iter = 1,
                       design = opdf[1:10, c(pars, y.name)] 
  )
  
  # with y specified as caharacter --> throws an error because y is chr but not in seen.points
  expect_error(
    p1.chr = PredictorAf$new(model = res$models[[1]],
                             data = res$seen.points[[1]], 
                             res.mbo = res,
                             iter = 1,
                             y = "ei",
                             design = opdf[1:10, c(pars, y.name)] 
    )
  )
  # with y comtained in the data
  ei = acq.fun(
    points = res$seen.points[[1]],
    models = list(res$models[[1]]),
    control = res$control,
    par.set = res$opt.path$par.set,
    designs = list(opdf[1:10, c(pars, y.name)]),
    iter = 1,
    progress = NULL,
    attributes = FALSE
  )
  p1.data = PredictorAf$new(model = res$models[[1]],
                            data = cbind(res$seen.points[[1]], ei = ei),
                            res.mbo = res,
                            iter = 1, 
                            design = opdf[1:10, c(pars, y.name)])
  
  # with y as vector
  p1.vec = PredictorAf$new(model = res$models[[1]],
                           data = res$seen.points[[1]],
                           res.mbo = res,
                           iter = 1, 
                           y = ei,
                           design = opdf[1:10, c(pars, y.name)])
  
  expect_equivalent(
    p1$predict(pp.true[1, ]),
    opdf[11, infill]
  )
  expect_equivalent(
    p1.data$predict(pp.true[1, ]),
    opdf[11, infill]
  )
  expect_equivalent(
    p1.vec$predict(pp.true[1, ]),
    opdf[11, infill]
  )
  
})

test_that("extracts data (seen.points) correctly from the mbo object", {
  ctrl = makeMBOControl(y.name = "target", store.model.at = 1:3)
  ctrl = setMBOControlTermination(ctrl, iters = 2)
  ctrl = setMBOControlInfill(ctrl,
                             opt = "focussearchSavepts",
                             opt.focussearch.maxit = 20,
                             opt.focussearch.points = 50,
                             crit = makeMBOInfillCritEI()
  )
  res = mbo(objfun,
            design = initial.data,
            control = ctrl,
            show.info = FALSE
  )
  opdf = as.data.frame(res$opt.path)
  pars = names(res$x)
  y.name = res$control$y.name
  infill = res$control$infill.crit$id
  pp.true = opdf[which(opdf$dob > 0), c(pars, y.name)]
  
  #create the PredictorAfObject in iter 1
  p1 = PredictorAf$new(model = res$models[[1]],
                       res.mbo = res,
                       iter = 1,
                       design = opdf[1:10, c(pars, y.name)] 
  )

  #create the PredictorAf Object in iter 2
  p2 = PredictorAf$new(model = res$models[[2]],
                       res.mbo = res,
                       iter = 2,
                       design = opdf[1:11, c(pars, y.name)] 
  )
  expect_equal(p1$data$X, data.table::data.table(res$seen.points[[1]]))
  expect_equal(p2$data$X, data.table::data.table(res$seen.points[[2]]))
  
})

test_that("PredictorAf gives error if data or design contain NA", {
  ctrl = makeMBOControl(y.name = "target", store.model.at = 1:3)
  ctrl = setMBOControlTermination(ctrl, iters = 2)
  ctrl = setMBOControlInfill(ctrl,
                             opt = "focussearchSavepts",
                             opt.focussearch.maxit = 20,
                             opt.focussearch.points = 50,
                             crit = makeMBOInfillCritEI()
  )
  res = mbo(objfun,
            design = initial.data,
            control = ctrl,
            show.info = FALSE
  )
  opdf = as.data.frame(res$opt.path)
  pars = names(res$x)
  y.name = res$control$y.name
  infill = res$control$infill.crit$id
  seen = res$seen.point[[1]]
  # data contains NA
  dat = rbind(seen, NA)
  expect_error(
    p.dat = PredictorAf$new(
      model = res$models[[1]],
      data = dat,
      res.mbo = res,
      iter = 1,
      design = opdf[1:10, c(pars, y.name)]
    )
  )
  
  # design contains NA
  des = rbind(
    opdf[1:9, c(pars, y.name)],
    NA
  )
  expect_error(
    p.des = PredictorAf$new(
      model = res$models[[1]],
      data = seen,
      res.mbo = res,
      iter = 1,
      design = des
    )
  )
})


test_that("PredictorAf gives error if design is not correctly specified",{
  ctrl = makeMBOControl(y.name = "target", store.model.at = 1:3)
  ctrl = setMBOControlTermination(ctrl, iters = 2)
  ctrl = setMBOControlInfill(ctrl,
                             opt = "focussearchSavepts",
                             opt.focussearch.maxit = 20,
                             opt.focussearch.points = 50,
                             crit = makeMBOInfillCritEI()
  )
  res = mbo(objfun,
            design = initial.data,
            control = ctrl,
            show.info = FALSE
  )
  opdf = as.data.frame(res$opt.path)
  pars = names(res$x)
  y.name = res$control$y.name
  infill = res$control$infill.crit$id
  pp.true = opdf[which(opdf$dob > 0), c(pars, y.name)]
  
  # missing y.name column if the design
  expect_error(p = PredictorAf$new(
    model = res$models[[1]],
    data = res$seen.points[[1]],
    res.mbo = res,
    iter = 1,
    design = opdf[1:10, pars]
  ))
  
  # wrong numbers of rows
  expect_error(p.ei = PredictorAf$new(
    model = res$models[[2]],
    data = res$seen.points[[2]],
    res.mbo = res,
    iter = 2,
    design = opdf[1:10, c(pars, y.name)]
  ))
  
})
