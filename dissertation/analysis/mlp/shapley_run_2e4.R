################################################################################
######################### Performance of sample size   #########################
################################################################################
#- to the best of our knowledge there is no way to find an optimal sample size for the 
#  Shapley Approximation (e.g. see Molnar iml book)
#- we try to find the optimal sample size in a greedy way using the efficiency gap
#  of the approximation (see checkSampleSize for more details). 
#  We expect that with increasing sample.size the approximation gets better as the 
#  efficiency gap becomes smaller.
#- We want to use the same sample size for each iteration, hence we find the best sample size
#  for best predicted in the process and use it for all proposals


# mbo results
mbo = readRDS(paste0(getwd(),"/analysis/mlp_phoneme/mbo_phoneme_lambda_1.rds"))
#find the iter of the best y (it could be that the best y is the initial design)
opdf = as.data.frame(mbo$opt.path)
# best predicted proposal 
sm113 = mbo$models$`113`
des = opdf[which(opdf$dob > 0), 1:7] # select only the proposals
des$mean = predict(sm113, newdata = des)$data$response
best.props = which.min(des$mean) # 95L --> used as iter.interest in the benchmark experiment


### RESULTS ###
# comparing the results of the same explanation with different sample size
size.1e2 = ShapleyMBO(res.mbo = mbo, iter.interest = 95, sample.size = 100, contribution = TRUE)
check.1e2 = checkSampleSize(size.1e2)
size.1e3 = ShapleyMBO(res.mbo = mbo, iter.interest = 95, sample.size = 1000, contribution = TRUE)
check.1e3 = checkSampleSize(size.1e3)
size.5e3 = ShapleyMBO(res.mbo = mbo, iter.interest = 95, sample.size = 5000, contribution = TRUE)
check.5e3 = checkSampleSize(size.5e3)
size.1e4 = ShapleyMBO(res.mbo = mbo, iter.interest = 95, sample.size = 10000, contribution = TRUE)
check.1e4 = checkSampleSize(size.1e4)
size.15e3 = ShapleyMBO(res.mbo = mbo, iter.interest = 95, sample.size = 15000, contribution = TRUE)
check.15e3 = checkSampleSize(size.15e3)
size.2e4 = ShapleyMBO(res.mbo = mbo, iter.interest = 95, sample.size = 20000, contribution = TRUE)
check.2e4 = checkSampleSize(size.2e4)
#--> sample size 2e4 is the only one for which each efficiency gap (cb, mean, se) 
#    is small enough. We therefore use 2e4 to explain both iter 95 and desirability path

### TIME ###
# comparing the execution time of the same explanation with different sample size
time.best.mean = rbenchmark::benchmark(
  size.1e2 = ShapleyMBO(res.mbo = mbo, iter.interest = 95, sample.size = 100, contribution = TRUE),
  size.1e3 = ShapleyMBO(res.mbo = mbo, iter.interest = 95, sample.size = 1000, contribution = TRUE),
  size.5e3 = ShapleyMBO(res.mbo = mbo, iter.interest = 95, sample.size = 5000, contribution = TRUE),
  size.1e4 = ShapleyMBO(res.mbo = mbo, iter.interest = 95, sample.size = 10000, contribution = TRUE),
  size.15e3 = ShapleyMBO(res.mbo = mbo, iter.interest = 95, sample.size = 15000, contribution = TRUE),
  size.2e4 = ShapleyMBO(res.mbo = mbo, iter.interest = 95, sample.size = 20000, contribution = TRUE),
  replications = 1
)
# storing the results
perf.shapley.mlp = list(
  "1e2" = list(shapley = size.1e2, sample.size = check.1e2),
  "1e3" = list(shapley = size.1e3, sample.size = check.1e3),
  "5e3" = list(shapley = size.5e3, sample.size = check.5e3),
  "1e4" = list(shapley = size.1e4, sample.size = check.1e4),
  "15e3" = list(shapley = size.15e3, sample.size = check.15e3),
  "2e4" = list(shapley = size.2e4, sample.size = check.2e4),
  "time" = time.best.mean
)
saveRDS(perf.shapley.mlp, "perf_shapley_mlp.rds")
#- 2e4 is as expected the most computationally expensive (ca. 90 times slower
#  than the default), yet overall computation time for the single explanation
#  is affordable (less than 5 min)

###########################################################
####   Application Example - MLP Phoneme        ###########
###########################################################
#- after benchmarking different sample sizes we finally explain the mbo path 
#- here we apply ShapleyMBO_mclapply (parallel version of ShapleyMBO) on the mbo 
#  results of the phoneme data set
#- sample.size = 20000
#- iter.interest = NULL (all iterations), seed = 1 (default)
start = Sys.time()
shapley = ShapleyMBO_mclapply(mbo, sample.size = 20000, contribution = TRUE, no.cores = 4)
end = Sys.time()
res = list("shapley" = shapley, "time.required" = end - start)
saveRDS(res, "shapley_phoneme_2e4.rds")
