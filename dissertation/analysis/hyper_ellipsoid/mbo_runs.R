########################################################
################## Experiments  ########################
########################################################
#- we run 30 optimization problems for each lambda (1, 10)
#- we set the seed once at the beginning of each lambda experiment as suggested by
#  Morris et al in "Using simulation studies to evaluate statistical methods"
#- the authors also suggest to store rng state before and after the single experiments.
#  We follow the suggestion and then store the rng.states in a RDS file
#- for different lambdas we set the same initial starting seed in order to 
#  have identical initial BO stochastic settings (initial design and obj.fun)
#- paths of the same run differ in general because of the cb lambda, sm and infill optimizer
#- as termination condition we use max.evals 80 (20d = 20*4) as used by Picheny et al. 
#  "A benchmark of kriging-based infill criteria for noisy optimization" (2013)
#- the authors also suggest to set the variance of the noise as a percentage of the std.dev. of the noise free
#  function. To have a moderate noise they suggest 5% --> sd(eps) = 0.05 * sd(noise.free.fun)
#- other arguments of createMBOrunHypEll4d_noisy are left to default values (opt = "focussearch", 
#  id.dim = 4, percent.noise = 0.05, seed = NULL)

# lambda1
rng.lambda1 = list()
set.seed(1)
for (i in 1:30) { # repeat 30 times
  oldseed = .GlobalEnv$.Random.seed
  createMBOrunHypEll4d_noisy(lambda = 1, run = i)
  newseed = .GlobalEnv$.Random.seed
  rng.lambda1[[i]] = list(oldseed = oldseed, newseed = newseed)
}
saveRDS(rng.lambda1, "rng_lambda1.rds")

# lambda10
rng.lambda10 = list()
set.seed(1)
for (i in 1:30) {
  oldseed = .GlobalEnv$.Random.seed
  createMBOrunHypEll4d_noisy(lambda = 10, run = i)
  newseed = .GlobalEnv$.Random.seed
  rng.lambda10[[i]] = list(oldseed = oldseed, newseed = newseed)
}
saveRDS(rng.lambda10, "rng_lambda10.rds")

# store the opdf in a list and then bind the elements of the list a single data frame
#lambda1
opdf_lambda_1 = list()
for (i in 1:30) { #loop over all bo runs with lambda1
  rds.file = paste0(getwd(),sprintf("/analysis/hyper_ellipsoid/mbo_runs/lambda_1/lambda_1_run_%i.rds", i))
  mbo = readRDS(rds.file)
  opdf = as.data.frame(mbo$opt.path)
  opdf_lambda_1[[i]] = opdf
}
opdf_lambda_1 = dplyr::bind_rows(opdf_lambda_1, .id = "run")
saveRDS(opdf_lambda_1, "opdf_lambda_1.rds")

#lambda10
opdf_lambda_10 = list()
for (i in 1:30) {
  rds.file = paste0(getwd(),sprintf("/analysis/hyper_ellipsoid/mbo_runs/lambda_10/lambda_10_run_%i.rds", i))
  mbo = readRDS(rds.file)
  opdf = as.data.frame(mbo$opt.path)
  opdf_lambda_10[[i]] = opdf
}
opdf_lambda_10 = dplyr::bind_rows(opdf_lambda_10, .id = "run")
saveRDS(opdf_lambda_10, "opdf_lambda_10.rds")

# bind the results together
opdf_lambda_all = dplyr::bind_rows(opdf_lambda_1, opdf_lambda_10, .id = "lambda") %>% 
  dplyr::mutate(lambda = factor(lambda, levels = c("1", "2"), labels = c("1", "10")))
saveRDS(opdf_lambda_all, "opdf_lambda_all.rds")
