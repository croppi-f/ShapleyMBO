################################################################################
##################   Apply ShapleyMBO on the MBO runs    #######################
################################################################################
# Shapley Settings:
# - res.mbo: every mbo results individually
# - iter.interest = NULL (default) --> all the 64 iterations
# - sample.size = 1000
# - contribution = TRUE --> apply the Lin.Axiom
# - seed = 1 (default)

no.cores = 4 # can be run also locally

# lambda 1
start.l1 = Sys.time()
shapley1e3_lambda_1 = parallel::mclapply(1:30, function(x) { # for all 30 BO runs
  # read the file
  rds.file = paste0(getwd(),sprintf("/dissertation/analysis/hyper_ellipsoid/mbo_runs_data/lambda_1/lambda_1_run_%i.rds", x))
  mbo = readRDS(rds.file)
  # apply ShapleyMBO (not the parallel version thereof)
  ShapleyMBO(mbo, sample.size = 1000, contribution = TRUE) # use default seed = 1
  },
  mc.cores = no.cores
)
end.l1 = Sys.time()
time.l1 = end.l1 - start.l1
shapley1e3_lambda_1 = dplyr::bind_rows(shapley1e3_lambda_1, .id = "run") # bind the results
saveRDS(shapley1e3_lambda_1, "shapley1e3_lambda_1.rds")

# lambda 10
start.l10 = Sys.time()
shapley1e3_lambda_10 = parallel::mclapply(1:30, function(x) {
  rds.file = paste0(getwd(),sprintf("/dissertation/analysis/hyper_ellipsoid/mbo_runs_data/lambda_10/lambda_10_run_%i.rds", x))
  mbo = readRDS(rds.file)
  ShapleyMBO(mbo, sample.size = 1000, contribution = TRUE) # use default seed = 1
},
mc.cores = no.cores
)
end.l10 = Sys.time()
time.l10 = end.l10 - start.l10
shapley1e3_lambda_10 = dplyr::bind_rows(shapley1e3_lambda_10, .id = "run") # bind the results
saveRDS(shapley1e3_lambda_10, "shapley1e3_lambda_10.rds")

# bind the Shapley results together in one sigle df
shapley1e3 = dplyr::bind_rows(shapley1e3_lambda_1, shapley1e3_lambda_10, .id = "lambda") %>% 
  dplyr::mutate(lambda = factor(lambda, levels = c("1", "2"), labels = c("1", "10")))
saveRDS(shapley1e3, "shapley1e3.rds")
