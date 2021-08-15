source("R/libraries.R")
################################################################################
##################   Apply ShapleyMBO on the MBO runs    #######################
################################################################################
# Shapley Settings:
# - res.mbo --> use all the 60 BO runs 
# - iter.interest = NULL (default) --> all the 64 iterations
# - sample.size = 1000
# - contribution = TRUE --> apply the Lin.Axiom
# - seed = 1 (default)

no.cores = 4 # run only on the LRZ Server with 10 cores

# lambda 1
start.l1 = Sys.time()
shapley1e3_lambda_1 = parallel::mclapply(1:30, function(x) { # for all 30 BO runs
  # read the file
  rds.file = paste0(getwd(),sprintf("/analysis/hyper_ellipsoid/mbo_runs/lambda_1/lambda_1_run_%i.rds", x))
  mbo = readRDS(rds.file)
  # apply ShapleyMBO
  ShapleyMBO(mbo, sample.size = 1000, contribution = TRUE) # use default seed = 1
  },
  mc.cores = no.cores
)
end.l1 = Sys.time()
time.l1 = end.l1 - start.l1 #1.369043 hours
shapley1e3_lambda_1 = dplyr::bind_rows(shapley1e3_lambda_1, .id = "run")
saveRDS(shapley1e3_lambda_1, "shapley1e3_lambda_1.rds")

# lambda 10
start.l10 = Sys.time()
shapley1e3_lambda_10 = parallel::mclapply(1:30, function(x) {
  rds.file = paste0(getwd(),sprintf("/analysis/hyper_ellipsoid/mbo_runs/lambda_10/lambda_10_run_%i.rds", x))
  mbo = readRDS(rds.file)
  ShapleyMBO(mbo, sample.size = 1000, contribution = TRUE) # use default seed = 1
},
mc.cores = no.cores
)
end.l10 = Sys.time()
time.l10 = end.l10 - start.l10 #2.767766 hours
shapley1e3_lambda_10 = dplyr::bind_rows(shapley1e3_lambda_10, .id = "run")
saveRDS(shapley1e3_lambda_10, "shapley1e3_lambda_10.rds")


