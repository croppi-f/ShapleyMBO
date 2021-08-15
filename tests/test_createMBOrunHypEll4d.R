######################### createMBOrunHypEll4d ##############################
######################### TEST ##############################################
library(testthat)

# we change the last lines of createMBOrunHypEll4d to check if the method works correctly
# (remore storing of object & create a result object instead)
createMBOrunHypEll4d_testthat = function(lambda,
                                         run,
                                         opt = "focussearch", 
                                         term.cond = "max.evals",
                                         max.evals = 80, #20d
                                         iters = 50, # enough to show convergence
                                         id.dim = 4,
                                         percent.noise = 0.05,
                                         seed = NULL) {
  # set.seed as first step of the work flow
  if (!is.null(seed))
    set.seed(seed) 
  # approximate the sd of the function
  noise.free.fun = makeHyperEllipsoidFunction(dimensions = 4L)
  sample = generateDesign(n = 10000, par.set = getParamSet(noise.free.fun), fun = lhs::randomLHS)
  values = apply(sample, 1, noise.free.fun)
  sd = sd(values)
  # sd of the noise is 5% of the estimated sd of the function (moderate noise)
  noise.sd = percent.noise * sd
  
  # define noisy Objective Function
  obj.fun = makeSingleObjectiveFunction(
    name = "noisy 4d Hyper-Ellipsoid",
    id = "hyper_ellipsoid_4d_5%noise",
    description = "4d Hyper-Ellipsoid with artificially added gaussian noise.
    The sd of the noise is 5% of sd of the noise free HypEll, eps ~ N(0, 0.05 * sd(nf.fun))",
    fn = function(x, sd = noise.sd) { #see makeHyperEllipsoidFunction()
      n = length(x)
      eps = rnorm(1, 0, sd) # Gaussian noise
      sum(1:n * x^2) + eps # the fist term is taken from the source code of makeHyperEllipsoidFunction()
    }, 
    par.set = makeNumericParamSet(
      len = 4, id = "x", # 4 dimensional 
      lower = rep(-5.12, 4), upper = rep(5.12, 4),
      vector = TRUE),
    global.opt.params = rep(0, 4), 
    global.opt.value = 0
  )
  # get the Parameter Set
  ps = getParamSet(obj.fun)
  
  # generate Initial Design
  size = id.dim * getNumberOfParameters(obj.fun) #4d
  des = generateDesign(n = size, par.set = ps, fun = lhs::maximinLHS)
  
  # set the Control object
  if (term.cond == "max.evals") dob = max.evals - size
  if (term.cond == "iters") dob = iters
  ctrl = makeMBOControl(store.model.at = 1:(dob + 1))# store every surrogate model
  ctrl = setMBOControlInfill(ctrl, crit = makeMBOInfillCritCB(cb.lambda = lambda), opt = opt)
  if (term.cond == "max.evals") ctrl = setMBOControlTermination(ctrl, max.evals = max.evals)
  if (term.cond == "iters") ctrl = setMBOControlTermination(ctrl, iters = iters)
  
  # 7. run the optimization problem
  mbo = mbo(obj.fun, design = des, control = ctrl, show.info = FALSE)
  opdf = as.data.frame(mbo$opt.path)
  
  res = list(
    init.des = opdf %>% dplyr::select(1:5) %>% dplyr::slice(1:size), # initial design
    # we exclude colums including "time" (exec.time, train.time, propose.time), since they are not really relevant this tests
    opdf = opdf[-c(9, 12, 14)], 
    # we search for noise.sd within the mbo results
    noise.sd = environment(environment(mbo[["final.opt.state"]][["opt.problem"]][["fun"]])[["fn"]])[["noise.sd"]],
    # we store the sm names is order to control we store the sm correctly using dob = max.evals - iters
    names.sm = names(mbo$models)
  )
             
}

# Setting - seed is set outside the function
# lambda1all
lambda1 = list()
set.seed(1)
for (i in 1:2) {
  lambda1[[i]] = createMBOrunHypEll4d_testthat(lambda = 1, term.cond = "max.evals", max.evals = 20)
}

lambda1_rep2 = list()
set.seed(1)
for (i in 1:2) {
  lambda1_rep2[[i]] = createMBOrunHypEll4d_testthat(lambda = 1, term.cond = "max.evals", max.evals = 20)
}

# lambda10
lambda10 = list()
set.seed(1)
for (i in 1:2) {
  lambda10[[i]] = createMBOrunHypEll4d_testthat(lambda = 10, term.cond = "max.evals", max.evals = 20)
}

######################### TEST ##############################################
test_that("createMBOrunHypEll4d is reproducible", {
  
  expect_identical(lambda1, lambda1_rep2)
  
}) #passed

test_that("different lambda have same obj.fun and init.des", {
  
  # obj.fun - same sd of the noise
  # fist run
  expect_identical(lambda1[[1]][["noise.sd"]], lambda10[[1]][["noise.sd"]])
  # second run
  expect_identical(lambda1[[2]][["noise.sd"]], lambda10[[2]][["noise.sd"]])
  
  # initial design
  # fist run
  expect_identical(lambda1[[1]][["init.des"]], lambda10[[1]][["init.des"]])
  # second run
  expect_identical(lambda1[[2]][["init.des"]], lambda10[[2]][["init.des"]])
  
}) # passed

test_that("all necessary surrogate models are stored - dob works correctly", {
  # bo is run with max.evals = 20. Therefore we expect to have in total 4 (20-16) iterations
  # or dob. Since store.model.at = 1:(dob + 1) we expect to have stored all the sm plus the final one
  
  expect_identical(lambda1[[1]]$names.sm, as.character(1:5))
  expect_identical(lambda1[[2]]$names.sm, as.character(1:5))
  expect_identical(lambda10[[1]]$names.sm, as.character(1:5))
  expect_identical(lambda10[[2]]$names.sm, as.character(1:5))
}) #passed
