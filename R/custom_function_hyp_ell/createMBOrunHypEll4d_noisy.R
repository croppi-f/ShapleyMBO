###########################################################
#########   createMBOrunHypEll4d   ########################
###########################################################
#- the function creates and stores BO runs on the Hyper Ellipsoid 4d with noise
#- the objective function is noisy with Gaussian noise, the std.dev of the noise 
#  is set as % of the std.dev. of the noise free 4d HypEll
#- the std.dev of the noise free 4d HypEll is approximated using a random (lhs::randomLHS)
#  sample of 10000 instances from the input space 
#- id.dim ("initial design dimension") is used to generate the initial design. id.dim = 4 --> a design
#  of size 4d is generated, with d being the dimension of the input space
createMBOrunHypEll4d_noisy = function(lambda,
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
  noise.free.fun = smoof::makeHyperEllipsoidFunction(dimensions = 4L)
  sample = ParamHelpers::generateDesign(n = 10000, par.set = getParamSet(noise.free.fun), fun = lhs::randomLHS)
  values = apply(sample, 1, noise.free.fun)
  sd = sd(values)
  # sd of the noise is 5% of the estimated sd of the function (moderate noise)
  noise.sd = percent.noise * sd
  
  # define noisy Objective Function
  obj.fun = smoof::makeSingleObjectiveFunction(
    name = "noisy 4d Hyper-Ellipsoid",
    id = "hyper_ellipsoid_4d_5%noise",
    description = "4d Hyper-Ellipsoid with artificially added gaussian noise.
    The sd of the noise is 5% of sd of the noise free HypEll, eps ~ N(0, 0.05 * sd(nf.fun))",
    fn = function(x, sd = noise.sd) { #see makeHyperEllipsoidFunction() of the smoof package
      n = length(x)
      eps = rnorm(1, 0, sd) # Gaussian noise
      sum(1:n * x^2) + eps # the fist term is taken from the source code of makeHyperEllipsoidFunction()
    }, 
    par.set = makeNumericParamSet(
      len = 4, id = "x", # 4 dimensional 
      lower = rep(-5.12, 4), upper = rep(5.12, 4),
      vector = TRUE),
    global.opt.params = rep(0, 4),
    global.opt.value = 0,
    noisy = TRUE
  )
  #Parameter Set
  ps = smoof::getParamSet(obj.fun)
  
  #Initial Design
  size = id.dim * smoof::getNumberOfParameters(obj.fun) #4 * 4 = 16
  des = ParamHelpers::generateDesign(n = size, par.set = ps, fun = lhs::maximinLHS)
  
  #Control object
  if (term.cond == "max.evals") dob = max.evals - size
  if (term.cond == "iters") dob = iters
  ctrl = makeMBOControl(store.model.at = 1:(dob + 1))# store every surrogate model
  ctrl = setMBOControlInfill(ctrl, crit = makeMBOInfillCritCB(cb.lambda = lambda), opt = opt)
  if (term.cond == "max.evals") ctrl = setMBOControlTermination(ctrl, max.evals = max.evals)
  if (term.cond == "iters") ctrl = setMBOControlTermination(ctrl, iters = iters)
  
  # run the optimization problem
  res = mbo(obj.fun, design = des, control = ctrl, show.info = FALSE)
  # store the results
  store_path = sprintf("lambda_%i_run_%i.rds", lambda, run)
  saveRDS(res, store_path)
}
