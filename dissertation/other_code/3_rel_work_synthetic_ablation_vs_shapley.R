############# Ablation vs. Shapley Value ###############
library(mlr)
library(iml)

# specify the data
set.seed(1) # for reproducibility
n = 10000
generateY = function(X) {
  eps = rnorm(nrow(X), sd = 0.05)
  form =  ~ V1 + V2:V3 - 1
  mat = model.matrix(form, data = X)
  rowSums(mat) + eps
}

V1 = runif(n, 0, 1)
V2 = runif(n, 0, 1)
V3 = runif(n, 0, 1)
X = as.data.frame(cbind(V1, V2, V3))
X$y = generateY(X = X)

# rows = sort(sample(1:10000, 8000))
# diff = setdiff(1:10000, rows)
# train_set = X[rows, ]
# test_set = X[diff, ]

# define the learning task
task = makeRegrTask(data = X, target = "y")
# define the learner
lrn = makeLearner("regr.randomForest", ntree = 100)
# train the learner
mod = mlr::train(lrn, task)

######################################################  
############## Shapley Value        ##################
######################################################
# instance to explain, which correspond to the target configuration
target = as.data.frame(t(c(0,0,0)))
#target$pred = predict(mod, newdata = target)$data$response
pred.target = predict(mod, newdata = target)$data$response

samples = X
samples$y = NULL


# bind explicand and samples
X.shapley = rbind(target, samples, make.row.names = FALSE)

# create the Predictor
P = Predictor$new(
  model = mod, 
  data = X.shapley
)

# create the Shapley object
S = Shapley$new(
  predictor = P,
  x.interest = X.shapley[1,],
  sample.size = 1000
)

payout.shapley = S$y.hat.interest - S$y.hat.average
names(payout.shapley) = NULL

shapley = S$results[, 1:2]
shapley$rel.imp = shapley$phi / payout.shapley
# note that the sum of the relative importance do not sum up to 1. This problem
# does not change the results and is problably caused by the approximation error

######################################################
############## Ablation Analysis #####################
###################################################### 
#- this time we try to find a instance which has similar values for V2 and V3 and 
#  and the same time a similar prediction to the average prediction in order to have
#  similar payputs for both shapley and ablation
source = as.data.frame(t(c(0.5,0.5,0.5)))
pred.source = predict(mod, newdata = source)$data$response
# payout ablation
payout.ablation = pred.source - pred.target
################ round one #######################
# flip V1
r1.v1 = source
r1.v1[["V1"]] = target[["V1"]]
pred.r1.v1 = predict(mod, newdata = r1.v1)$data$response
gain.r1.v1 = pred.source - pred.r1.v1
# flip V2
r1.v2 = source
r1.v2[["V2"]] = target[["V2"]]
pred.r1.v2 = predict(mod, newdata = r1.v2)$data$response
gain.r1.v2 = pred.source - pred.r1.v2
# flip V3
r1.v3 = source
r1.v3[["V3"]] = target[["V3"]]
pred.r1.v3 = predict(mod, newdata = r1.v3)$data$response
gain.r1.v3 = pred.source - pred.r1.v3

# determine which param to flip
r1 = c(
  "V1" = gain.r1.v1,
  "V2" = gain.r1.v2,
  "V3" = gain.r1.v3
)
best.r1 = names(r1)[which.max(r1)] # "v1" --> flip V1

################ round two #######################
#- now the reference configuration becomes r1.v1, since we flipped
#  the value of V1 to its target value in the precedent round
#flip V2
r2.v2 = r1.v1
r2.v2[["V2"]] = target[["V2"]]
pred.r2.v2 = predict(mod, newdata = r2.v2)$data$response
gain.r2.v2 = pred.r1.v1 - pred.r2.v2
#flip V3
r2.v3 = r1.v1
r2.v3[["V3"]] = target[["V3"]]
pred.r2.v3 = predict(mod, newdata = r2.v3)$data$response
gain.r2.v3 = pred.r1.v1 - pred.r2.v3
# determine which param to flip
r2 = c(
  "V2" = gain.r2.v2,
  "V3" = gain.r2.v3
)
best.r2 = names(r2)[which.max(r2)] # "v2" --> flip V2

################ round three #####################
#- in the final round only V3 is left to flip
gain.r3.v3 = pred.r2.v2 - pred.target
####Putting up Ablation results together
marg.contr = c(
  r1[best.r1], 
  r2[best.r2],
  "V3" = gain.r3.v3
)
sum(marg.contr) == payout.ablation #TRUE

ablation = data.frame(
  feature = names(marg.contr),
  standalone = r1,
  mc = marg.contr,
  rel.imp = marg.contr / payout.ablation,
  row.names = as.character(1:3)
)

flipping.order = c(best.r1, best.r2, "V3")

####Putting up Shapley and Ablation together
final.res = dplyr::left_join(
  ablation[, c(1,2,4)],
  shapley[,c(1,3)],
  by = "feature",
  suffix = c(".ablation", ".shapley")
) 

 