############################################
###    Motivating Example             ######
############################################
# univariate
# set the seed for reproducibility 
set.seed(1)
obj.fun = makeCosineMixtureFunction(1)
obj.fun = convertToMinimization(obj.fun)
des = generateDesign(n = 4, par.set = getParamSet(obj.fun), fun = lhs::randomLHS)
control = makeMBOControl()
control = setMBOControlTermination(control, iters = 3)
control = setMBOControlInfill(control, crit = makeMBOInfillCritCB(cb.lambda = 1))
mbo = exampleRun(obj.fun, design = des, control = control, show.info = FALSE)
plot1 = plotExampleRun(
  mbo, 
  iters = 1,
  pause = TRUE,
  point.size = 5,
  gg.objects = list(
    theme_bw(),
    theme(text = element_text(size = 20),
          legend.position = "none"
    )
  )
)
plot2 = plotExampleRun(
  mbo, 
  iters = 2,
  pause = TRUE,
  point.size = 5,
  gg.objects = list(
    theme_bw(),
    theme(text = element_text(size = 20),
          legend.position = "none"
    )
  )
)
plot3 = plotExampleRun(
  mbo, 
  iters = 3,
  pause = TRUE,
  point.size = 5,
  gg.objects = list(
    theme_bw(),
    theme(text = element_text(size = 20),
          legend.position = "bottom"
    )
  )
)

# 4 dimensional
set.seed(1)
obj.fun = makeCosineMixtureFunction(3)
obj.fun = convertToMinimization(obj.fun)
des = generateDesign(n = 4, par.set = getParamSet(obj.fun), fun = lhs::randomLHS)
control = makeMBOControl()
control = setMBOControlTermination(control, iters = 3)
control = setMBOControlInfill(control, crit = makeMBOInfillCritCB(cb.lambda = 1))
mbo = mbo(obj.fun, design = des, control = control, show.info = FALSE)
opdf = as.data.frame(mbo$opt.path)
proposed = opdf[5:7, c(1:3,4,9,14,15)]
