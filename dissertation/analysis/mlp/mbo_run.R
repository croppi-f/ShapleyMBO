###########################################################
####   Application Example - MLP Phoneme        ###########
###########################################################
# seed for reproducibility
set.seed(1)

# Parameter Set
ps = makeParamSet(
  makeNumericParam("batch_size", lower = log(16, 2), upper = log(512, 2), trafo = function(x) round(2^x)), 
  makeNumericParam("max_dropout", lower = 0, upper = 1), 
  makeNumericParam("max_units", lower = log(64, 2), upper = log(1024, 2), trafo = function(x) round(2^x)), 
  makeIntegerParam("num_layers", lower = 1, upper = 5),
  makeNumericParam("learning_rate", lower = 0, upper = 0.01),
  makeNumericParam("momentum", lower = 0.1, upper = 1),
  makeNumericParam("weight_decay", lower = 0, upper = 0.1)
)

# Objective Function
# read the file, where the EPM is stored
surr = readRDS(paste0(getwd(),"/dissertation/analysis/mlp/mbo_run_data/surrogate.rds"))
# EPMs - used as a surrogate target function (surr.val) and to evaluate the configuration
# on the test data (surr.test)
surr.val = surr$result[[1]][["model_val_balanced_acc"]][[1]]
surr.test = surr$result[[1]][["model_test_balanced_acc"]][[1]]

obj.fun = makeSingleObjectiveFunction(
  name = "mlp.surr.tuning",
  fn = function(x) {
    x = as.data.frame(x)
    y = predict(surr.val, newdata = x)$data$response
    attr(y, "extras") = list(test_performance = predict(surr.test, newdata = x)$data$response)
    return(1 - y / 100)
  },
  par.set = ps,
  noisy = TRUE,
  has.simple.signature = FALSE,
  minimize = TRUE
)

# Initial Design
size = 4 * getNumberOfParameters(obj.fun) #4d
des = generateDesign(n = size, par.set = ps, fun = lhs::maximinLHS)

# Control Flow
max.evals = 20 * getNumberOfParameters(obj.fun) #20d
dob = max.evals - size
ctrl = makeMBOControl(store.model.at = 1:(dob + 1))# store every surrogate model
ctrl = setMBOControlInfill(ctrl, crit = makeMBOInfillCritCB(cb.lambda = 1), opt = "focussearch") #cb with lambda = 1
ctrl = setMBOControlTermination(ctrl, max.evals = max.evals)

# run the optimization problem
res = mbo(obj.fun, design = des, control = ctrl, show.info = FALSE)
saveRDS(res, "mbo_phoneme_lambda_1.rds")