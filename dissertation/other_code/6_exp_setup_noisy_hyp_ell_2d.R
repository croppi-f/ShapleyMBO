library(smoof)
###############################
############# Bivariate Noisy Hyper Ellipsoid
###############################
#- here  we want to show how a possible noisy Hyp Ell can look like (figure 2 in thesis)
set.seed(1)
# approximate the sd of the function
noise.free.fun = makeHyperEllipsoidFunction(dimensions = 2L)
sample = generateDesign(n = 10000, par.set = getParamSet(noise.free.fun), fun = lhs::randomLHS)
values = apply(sample, 1, noise.free.fun)
sd = sd(values)
# sd of the noise is 5% of the estimated sd of the function (moderate noise)
noise.sd = 0.05 * sd

noisyHypEll2d = makeSingleObjectiveFunction(
  name = "",
  id = "hyper_ellipsoid_2d_5%noise",
  description = "2d Hyper-Ellipsoid with artificially added gaussian noise.
    The sd of the noise is 5% of the noise free function sd. eps ~ N(0, 0.05 * sd(fun))",
  fn = function(x, sd = noise.sd) {
    n = length(x)
    eps = rnorm(1, 0, sd) # Gaussian noise
    sum(1:n * x^2) + eps 
  }, 
  par.set = makeNumericParamSet(
    len = 2, id = "x",
    lower = rep(-5.12, 2), upper = rep(5.12, 2),
    vector = TRUE),
  global.opt.params = rep(0, 2), 
  global.opt.value = 0
)

plot(noisyHypEll2d, render.contours = TRUE, render.levels = TRUE)

#- this example is used for the expectation on the mean contribution (footnote in section 6.1)
# what is the 4d noise free hyp ell change if we replace a similar parameter with the average distance value from its optimum?
# we omit using noise for  the example
#- note that higher pars more distant from their optimum & nevertheless have a bigger impact.
hypEll = makeHyperEllipsoidFunction(dimensions = 4L)
y = hypEll(c(0.05, 0.15, 0.25, 0.35))
y.replaced= c(
  y.no1 = hypEll(c(2.56, 0.15, 0.25, 0.35)),
  y.no2 = hypEll(c(0.05, 2.56, 0.25, 0.35)),
  y.no3 = hypEll(c(0.05, 0.15, 2.56, 0.35)),
  y.no4 = hypEll(c(0.05, 0.15, 0.25, 2.56))
)
y.change = round(y - y.replaced, 2)
