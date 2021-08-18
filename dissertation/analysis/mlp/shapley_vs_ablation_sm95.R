############################################
###    Comparison with Ablation       ######
############################################
#- here we compare the mean contribution of the best iteration
#  found in the process (iter95) with the contributions of the 
#  same instance that would result using Ablation Analysis
#- to make results the most comparable we use as source configuration the instance with
#  prediction most similar to the average prediction in order to have a  
#  similar payout
#- as EPM for the ablation analysis we use the bo surrogate model 
#  of iter 95
#- NOTE: since we do not have found any R package yet that allow a flexible used of the 
#  ablation analysis algorithm we conduct the analysis ourselves. 
#  (e.g. teh frace package requires to tune the model with the frace package to conduct ablation analysis)
#  For that we followed the procedure give in "Analysing differences between algorithm configurations through ablation" 
#  by Fawcett and Hoos.

# set seed for reproducibility
set.seed(1)

mbo = readRDS(paste0(getwd(),"/dissertation/analysis/mlp/mbo_run_data/mbo_phoneme_lambda_1.rds"))
opdf = as.data.frame(mbo$opt.path)
sm = mbo$models$`95` # extract the surrogate model fot the ablation analysis
shapley = readRDS(paste0(getwd(), "/dissertation/analysis/mlp/shapley_run_data/shapley_phoneme_2e4.rds"))
# take only the results and not the time required
shapley = shapley[[1]]
################################
######## Shapley Results #######
################################
s95 = shapley[which(shapley$iter == "95"), ]
payout.mean = unique(s95$pred.interest_mean - s95$pred.average_mean)
s95$importance = s95$phi_mean / payout.mean # relative importance 
shapley.imp = s95[, c(2,18)]
rownames(shapley.imp) = NULL
################################
######## Ablation Analysis #####
################################
#- we measure the marginal contribution of parameters in terms of "gain", that is
#  we start by the source configuration and sequentially flip all pars until
#  we get to the target configuration.
#- in each round we flip / change the value of that parameter associated with the
#  biggest gain. gain > 0 indicates that modified configuration performs better
#  than the unmodified. Since we the objective function in MBO is minimized the lower the prediction 
#  of modified configurations, teh bigger the gain & the bigger and better their marginal contributions

#find source configuration
avg.pred = unique(s95$pred.average_mean) # reference value
# sample the obs from the space
samples = ParamHelpers::generateDesign(n = 7000, par.set = mbo$opt.path$par.set, fun = lhs::randomLHS)
samples$mean = predict(sm, newdata = samples)$data$response # predict 
samples$delta = abs(samples$mean - avg.pred) #the distance to the sv reference value
source = samples[which.min(samples$delta), 1:7] # select the instance which is closest
rownames(source) = NULL
pred.source = predict(sm, newdata = source)$data$response # source prediction
# target configuration
target = opdf[which(opdf$dob == 95), 1:7] # same as the instance explained by ShapleyMBO
rownames(target) = NULL
pred.target = predict(sm, newdata = target)$data$response # target prediction
explain = pred.source - pred.target # quantity to be explained (> 0 as expected, means that target has better 
# "performance" than source, recall that y is minimized)

#### START THE ABLATION ANALYSIS
################ round one #######################
# batch_size, max_dropout, max_units, num_layers, learning_rate, momentum, weight_decay 
# are the pars to flip
# flip batch_size
r1.bs = source
r1.bs[["batch_size"]] = target[["batch_size"]]
pred.r1.bs = predict(sm, newdata = r1.bs)$data$response
gain.r1.bs = pred.source - pred.r1.bs
# flip max_dropout
r1.md = source
r1.md[["max_dropout"]] = target[["max_dropout"]]
pred.r1.md = predict(sm, newdata = r1.md)$data$response
gain.r1.md = pred.source - pred.r1.md
# flip max_units
r1.mu = source
r1.mu[["max_units"]] = target[["max_units"]]
pred.r1.mu = predict(sm, newdata = r1.mu)$data$response
gain.r1.mu = pred.source - pred.r1.mu
# flip num_layers
r1.nl = source
r1.nl[["num_layers"]] = target[["num_layers"]]
pred.r1.nl = predict(sm, newdata = r1.nl)$data$response
gain.r1.nl = pred.source - pred.r1.nl
# flip learning_rate
r1.lr = source
r1.lr[["learning_rate"]] = target[["learning_rate"]]
pred.r1.lr = predict(sm, newdata = r1.lr)$data$response
gain.r1.lr = pred.source - pred.r1.lr
# flip momentum
r1.mom = source
r1.mom[["momentum"]] = target[["momentum"]]
pred.r1.mom = predict(sm, newdata = r1.mom)$data$response
gain.r1.mom = pred.source - pred.r1.mom
# flip weight_decay
r1.wd = source
r1.wd[["weight_decay"]] = target[["weight_decay"]]
pred.r1.wd = predict(sm, newdata = r1.wd)$data$response
gain.r1.wd = pred.source - pred.r1.wd

# determine which param to flip definitely
r1 = c(
  "bs" = gain.r1.bs, 
  "md" = gain.r1.md, 
  "mu" = gain.r1.mu, 
  "nl" = gain.r1.nl, 
  "lr" = gain.r1.lr, 
  "mom" = gain.r1.mom, 
  "wd" = gain.r1.wd
)
best.r1 = names(r1)[which.max(r1)] # "nl" --> flip num_layers

################ round two   #######################
#- now the reference configuration becomes r1.nl, since we flipped
#  the value of num_layers to its target value in the precedent round
# flip max_dropout
r2.md = r1.nl
r2.md[["max_dropout"]] = target[["max_dropout"]]
pred.r2.md = predict(sm, newdata = r2.md)$data$response
gain.r2.md = pred.r1.nl - pred.r2.md
# flip max_units
r2.mu = r1.nl
r2.mu[["max_units"]] = target[["max_units"]]
pred.r2.mu = predict(sm, newdata = r2.mu)$data$response
gain.r2.mu = pred.r1.nl - pred.r2.mu
# flip batch_size
r2.bs = r1.nl
r2.bs[["batch_size"]] = target[["batch_size"]]
pred.r2.bs = predict(sm, newdata = r2.bs)$data$response
gain.r2.bs = pred.r1.nl - pred.r2.bs
# flip learning_rate
r2.lr = r1.nl
r2.lr[["learning_rate"]] = target[["learning_rate"]]
pred.r2.lr = predict(sm, newdata = r2.lr)$data$response
gain.r2.lr = pred.r1.nl - pred.r2.lr
# flip momentum
r2.mom = r1.nl
r2.mom[["momentum"]] = target[["momentum"]]
pred.r2.mom = predict(sm, newdata = r2.mom)$data$response
gain.r2.mom = pred.r1.nl - pred.r2.mom
# flip weight_decay
r2.wd = r1.nl
r2.wd[["weight_decay"]] = target[["weight_decay"]]
pred.r2.wd = predict(sm, newdata = r2.wd)$data$response
gain.r2.wd = pred.r1.nl - pred.r2.wd
# determine which param to flip definitely (num layers is not included!)
r2 = c(
  "md" = gain.r2.md, 
  "mu" = gain.r2.mu, 
  "bs" = gain.r2.bs, 
  "lr" = gain.r2.lr, 
  "mom" = gain.r2.mom, 
  "wd" = gain.r2.wd
)
best.r2 = names(r2)[which.max(r2)] # "wd" --> flip weight_decay

################ round three   #######################
#- now the reference configuration becomes r2.wd, since we also flipped
#  the value of weight_decay to its target value in r2. r2wd has wd and nl values flipped to their target and the
#  rest still to their source value
# flip max_dropout
r3.md = r2.wd
r3.md[["max_dropout"]] = target[["max_dropout"]]
pred.r3.md = predict(sm, newdata = r3.md)$data$response
gain.r3.md = pred.r2.wd - pred.r3.md
# flip max_units
r3.mu = r2.wd
r3.mu[["max_units"]] = target[["max_units"]]
pred.r3.mu = predict(sm, newdata = r3.mu)$data$response
gain.r3.mu = pred.r2.wd - pred.r3.mu
# flip batch_size
r3.bs= r2.wd
r3.bs[["batch_size"]] = target[["batch_size"]]
pred.r3.bs = predict(sm, newdata = r3.bs)$data$response
gain.r3.bs = pred.r2.wd - pred.r3.bs
# flip learning_rate
r3.lr = r2.wd
r3.lr[["learning_rate"]] = target[["learning_rate"]]
pred.r3.lr = predict(sm, newdata = r3.lr)$data$response
gain.r3.lr = pred.r2.wd - pred.r3.lr
# flip momentum
r3.mom = r2.wd
r3.mom[["momentum"]] = target[["momentum"]]
pred.r3.mom = predict(sm, newdata = r3.mom)$data$response
gain.r3.mom = pred.r2.wd - pred.r3.mom
# determine which param to flip definitely
r3 = c(
  "md" = gain.r3.md, 
  "mu" = gain.r3.mu, 
  "bs" = gain.r3.bs, 
  "lr" = gain.r3.lr, 
  "mom" = gain.r3.mom
)
best.r3 = names(r3)[which.max(r3)] # "bs" --> flip batch_size

################ round four   #######################
#- now the reference configuration becomes r3.bs, since we flipped
#  the value of batch_size to its target value in r3
# flip max_dropout
r4.md = r3.bs
r4.md[["max_dropout"]] = target[["max_dropout"]]
pred.r4.md = predict(sm, newdata = r4.md)$data$response
gain.r4.md = pred.r3.bs - pred.r4.md
# flip max_units
r4.mu = r3.bs
r4.mu[["max_units"]] = target[["max_units"]]
pred.r4.mu = predict(sm, newdata = r4.mu)$data$response
gain.r4.mu = pred.r3.bs - pred.r4.mu
# flip learning_rate
r4.lr = r3.bs
r4.lr[["learning_rate"]] = target[["learning_rate"]]
pred.r4.lr = predict(sm, newdata = r4.lr)$data$response
gain.r4.lr = pred.r3.bs - pred.r4.lr
# flip momentum
r4.mom = r3.bs
r4.mom[["momentum"]] = target[["momentum"]]
pred.r4.mom = predict(sm, newdata = r4.mom)$data$response
gain.r4.mom = pred.r3.bs - pred.r4.mom
# determine which param to flip definitely
r4 = c(
  "md" = gain.r4.md, 
  "mu" = gain.r4.mu, 
  "lr" = gain.r4.lr, 
  "mom" = gain.r4.mom
)
best.r4 = names(r4)[which.max(r4)] # "mu" --> flip max_units

################ round five   #######################
#- now the reference configuration becomes r4.mu, since we also flipped
#  the value of max_units to its target value in r4
# flip max_dropout
r5.md = r4.mu
r5.md[["max_dropout"]] = target[["max_dropout"]]
pred.r5.md = predict(sm, newdata = r5.md)$data$response
gain.r5.md = pred.r4.mu - pred.r5.md
# flip learning_rate
r5.lr = r4.mu
r5.lr[["learning_rate"]] = target[["learning_rate"]]
pred.r5.lr = predict(sm, newdata = r5.lr)$data$response
gain.r5.lr = pred.r4.mu - pred.r5.lr
# flip momentum
r5.mom = r4.mu
r5.mom[["momentum"]] = target[["momentum"]]
pred.r5.mom = predict(sm, newdata = r5.mom)$data$response
gain.r5.mom = pred.r4.mu - pred.r5.mom
# determine which param to flip
r5 = c(
  "md" = gain.r5.md, 
  "lr" = gain.r5.lr, 
  "mom" = gain.r5.mom
)
best.r5 = names(r5)[which.max(r5)] # "mom" --> flip momentum

################ round six   #######################
#- now the reference configuration becomes r5.mom, since we also flipped
#  the value of learning_rate to its target value in r5
# flip max_dropout
r6.md = r5.mom
r6.md[["max_dropout"]] = target[["max_dropout"]]
pred.r6.md = predict(sm, newdata = r6.md)$data$response
gain.r6.md = pred.r5.mom - pred.r6.md
# flip learning rate
r6.lr = r5.mom
r6.lr[["learning_rate"]] = target[["learning_rate"]]
pred.r6.lr = predict(sm, newdata = r6.lr)$data$response
gain.r6.lr = pred.r5.mom - pred.r6.lr
# determine which param to flip
r6 = c(
  "md" = gain.r6.md, 
  "lr" = gain.r6.lr
)
best.r6 = names(r6)[which.max(r6)] # "lr" --> flip learning_rate

################ round seven   #######################
#- in the final round only max_dropout is left to flip
gain.r7.md = pred.r6.lr - pred.target

################################################
############ Putting the results of AA together
################################################
# the final results
marg.contr= c(
  r1[best.r1], 
  r2[best.r2],
  r3[best.r3],
  r4[best.r4],
  r5[best.r5],
  r6[best.r6],
  "md" = gain.r7.md
)
sum(marg.contr) == explain #TRUE
ablation = data.frame(
  feature.id= names(marg.contr),
  mc = marg.contr,
  row.names = as.character(1:7)
)
flipping.order = c(best.r1, best.r2, best.r3, best.r4, best.r5, best.r6, "md")
ablation$feature = c("num_layers", "weight_decay", "batch_size",  "max_units", "momentum", "learning_rate", "max_dropout")
ablation = ablation[, c(3, 2)]
ablation$importance = ablation$mc / explain # relative imprtance to compare with the shapley results
ablation$round = 1:7 # blation round
ablation$cumsum = cumsum(ablation$mc) # cumulative sum of the marginal contributions

############################################################
############ Putting the results of Shapley and AA together
############################################################
# table
final.res = dplyr::left_join(
  shapley.imp,
  ablation[,c(1,3)],
  by = "feature",
  suffix = c(".shapley", ".ablation")
) 
#pivot to longer format to plot
final.res.long = final.res %>% tidyr::pivot_longer(
  cols = 2:3, 
  names_to = "method", 
  values_to = "importance", 
  names_prefix = "importance."
)
# plot settings
theme = theme_bw() + 
  theme(text = element_text(size = 15),
        legend.position = "bottom"
  )
order = forcats::fct_reorder(shapley$feature, shapley$phi_mean)
lev = levels(order)
# figure 11 i
plot.imp = ggplot(
  final.res.long,
  aes(x = forcats::fct_relevel(feature, lev), y = importance)
) +
  geom_col(fill = "#999999") +
  scale_x_discrete(labels = c("bs", "nl", "lr", "mu", "wd", "md", "mom")) +
  facet_grid(rows = vars(forcats::fct_relevel(method, c("shapley", "ablation")))) +
  labs(y = "relative.importance", x = "parameter") +
  theme 

############################################################
############ H-Stat      ##################################
############################################################
#- here we try to investigate if the surrogate models learned interactions effects,
#  in particular between LR and MU and other variables
#- the goal is to investigaet differences between the AA ans SV results
#- if so than we can effectively say than SV is a better choice  than Ablation

### IMPORTANT ###
#- consider whether to run or not the computations b/c they require a lot of time, especially the last one!
ilr = Interaction$new(P, grid.size = 100, feature = "learning_rate")
imu = Interaction$new(P, grid.size = 100, feature = "max_units") 
interactions = Interaction$new(P, grid.size = 100)

h.stat = list("all" = interactions, "learning_rate" = ilr, "max_units" = imu)
saveRDS(h.stat, "h_statistic.rds")

# single plots
h.all = plot(h.stat[[1]]) + theme + scale_y_discrete(labels = c("mom", "md", "wd", "mu", "bs", "nl", "lr")) + ylab("") # not used in the thesis
h.lr = plot(h.stat[[2]]) + theme + scale_y_discrete(labels = c("nl : lr", "mu : lr", "mom : lr", "md : lr", "wd : lr", "bs : lr")) + ylab("")
h.mu = plot(h.stat[[3]]) + theme + scale_y_discrete(labels = c("mom : mu", "lr : mu", "bs : mu", "md : mu", "wd : mu", "nl : mu")) + ylab("")
plot.h = (h.lr + h.mu) & xlim(0, 0.3) & xlab("interaction strength")

# figure 11 ii and iii
final.plot = plot.imp / plot.h + plot_annotation(tag_levels = "i")
