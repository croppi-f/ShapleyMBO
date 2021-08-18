library(tidyverse)
library(patchwork)
###########################################################
####   Application Example - MLP Phoneme  - Evaluation ####
###########################################################
# common layers for the plots and additional stuff
theme = theme_bw() + 
  theme(text = element_text(size = 15),
        legend.position = "bottom"
  )
hline = geom_hline(yintercept = 0, alpha = 0.5, linetype = "longdash")

################################################################################
###############         MBO results           ##################################
################################################################################
#loead the results
mbo = readRDS(paste0(getwd(),"/dissertation/analysis/mlp/mbo_run_data/mbo_phoneme_lambda_1.rds"))
opdf = as.data.frame(mbo$opt.path)
path = opdf[which(opdf$dob > 0), ] #proposals
rownames(path) = path$dob

# best predicted proposal 
sm113 = mbo$models$`113`
props = opdf[which(opdf$dob > 0), 1:7] # select only the proposals
props$mean = predict(sm113, newdata = props)$data$response
best.props = which.min(props$mean) # 95L
# best y of the process
best.y.id = which.min(opdf$y) #123L, which is iter 95L in path!
# best predicted including initial design
des = opdf[, 1:7]
des$mean = predict(sm113, newdata = des)$data$response
best.pred = which.min(des$mean) # 123L, which is iter 95L in path!
#--> proposal in iter 95 is also the best.y and best.predicted among all visited points! 
# same scale for y, mean and cb
lims = c(
  min(path$y, path$mean, path$cb),
  max(path$y, path$mean, path$cb)
)

# single prediction (error, mean, se, cb) plots, skip and go directly to summary.bo.mlp
{
  plot.y = ggplot(path, aes(x = dob, y = y)) +
    geom_line() +
    geom_point(aes(x = 95, y = opdf[which(opdf$dob == 95), "y"]), color = "#F8766D" ) +
    geom_point(aes(x = 108, y = opdf[which(opdf$dob == 108), "y"]), color = "#00BA38") +
    scale_y_continuous(limits = lims) +
    labs(x = "iter", y = "error") +
    theme
  plot.cb = ggplot(path, aes(x = dob, y = cb)) +
    geom_line() +
    geom_point(aes(x = 95, y = opdf[which(opdf$dob == 95), "cb"]), color = "#F8766D") +
    geom_point(aes(x = 108, y = opdf[which(opdf$dob == 108), "cb"]), color = "#00BA38") +
    scale_y_continuous(limits = lims) +
    labs(x = "iter") +
    theme
  plot.mean = ggplot(path, aes(x = dob, y = mean)) +
    geom_line() +
    geom_point(aes(x = 95, y = opdf[which(opdf$dob == 95), "mean"]), color = "#F8766D") +
    geom_point(aes(x = 108, y = opdf[which(opdf$dob == 108), "mean"]), color = "#00BA38") +
    scale_y_continuous(limits = lims) +
    labs(x = "iter") +
    theme
  plot.se = ggplot(path, aes(x = dob, y = se)) +
    geom_line() +
    geom_point(aes(x = 95, y = opdf[which(opdf$dob == 95), "se"]), color = "#F8766D") +
    geom_point(aes(x = 108, y = opdf[which(opdf$dob == 108), "se"]), color = "#00BA38") +
    labs(x = "iter") +
    theme
}
#figure 9
summary.bo.mlp = plot.y / plot.mean / plot.se / plot.cb

################################################################################
###############         Shapley results           ##############################
################################################################################
# colors: "#F8766D" "#53B400" "#A58AFF" "#00C094" "#FB61D7" "#C49A00" "#00B6EB"
# load the results
shapley = readRDS(paste0(getwd(), "/dissertation/analysis/mlp/shapley_run_data/shapley_phoneme_2e4.rds"))
# take only the results and not the time required
shapley = shapley[[1]]
# define the order of the feature
order = forcats::fct_reorder(shapley$feature, shapley$phi_cb)
lev = levels(order)

################################################################################
############# 1. Explain best predicted proposal ###############################
################################################################################
shapley95 = shapley[which(shapley$iter == "95"), ] # iter 95
# table 10 in thesis
table95 = t(shapley95[c(1, 4, 3, 5, 7, 2, 6),c(16, 8, 12)])
colnames(table95) = shapley95$feature[c(1, 4, 3, 5, 7, 2, 6)]

#figure 10
plot.best = plotShapleyMBO(
  shapley, infill.mbo = "cb", iter.interest = 95, lambda.mbo = 1, 
  sample.size = 20000, type = "bar", avg.phi = FALSE, ci = TRUE, 
  ci.alpha = 0.05, decomp = 3L
) & scale_y_continuous(limits = c(-0.02, 0.003)) # we change the scale b/c later on we also display iter 108
#- now we investigate if the results of the mean contribution would differ if we
#  Ablation Analysis is used instead of Shapley Value

# --> see .R script "shapley_vs_ablation_sm95"!

################################################################################
################# 2. Desirability paths     ####################################
################################################################################
# average contribution over the entire process
avg.phi = shapley %>% dplyr::select(feature, phi_mean, phi_se_scaled, phi_cb) %>%
  dplyr::group_by(feature) %>% 
  dplyr::summarise(avg.mean = mean(phi_mean), avg.se_scaled = mean(phi_se_scaled), avg.cb = mean(phi_cb))

# here we use the output of plotShapleyMBO as a base plot (all figure in one plot) but 
# then selected elements only to display plots on top of each other (see [[]])
base.plot = plotShapleyMBO(
  shapley, infill.mbo = "cb", iter.interest = NULL, lambda.mbo = 1, 
  sample.size = 20000, type = "line", avg.phi = FALSE, ci = TRUE, 
  ci.alpha = 0.05, decomp = 1L
)
cb.path = base.plot[[1]] & theme(legend.position = "none")
mean.path = base.plot[[2]][[1]] & theme(legend.position = "none")
se.path = base.plot[[2]][[2]]
#figure 12
plot.paths = cb.path / mean.path / se.path & 
  geom_vline(xintercept = 84L, linetype = "dotted") & 
  geom_vline(xintercept = 95L, linetype = "dotted") & 
  geom_vline(xintercept = 108L, linetype = "dotted")

# parallel plot of the configurations seen in the process (init des included)
# figure 13
par.plot = GGally::ggparcoord(opdf, columns = 1:7, groupColumn = 9, order = c(6,2,7,3,5,4,1))  + 
  scale_color_viridis_c(begin = 0, end = 1) +
  labs(x = "", y = "", color = "iter") +
  theme +
  scale_x_discrete(labels = c("mom", "md", "wd", "mu", "lr", "nl", "bs"))

################################################################################
########### 3. Learning rate final exploration         #########################
################################################################################
#- here we further investigate the final exploration rush of the learning rate dimension
# and take proposal of iter 108 as an example
shapley108 = dplyr::slice(shapley, which(iter == "108")) 
# highlighted parallel plot for iter 108 and iter 95 only:
#- we add a fake color variable in the data in order to highlight both iters with discrete colors
opdf$fake.color = c(rep("no.color", 122), "color1", rep("no.color", 12), "color2", rep("no.color", 4))
#figure 14 right
par.plot108 = GGally::ggparcoord(opdf, columns = 1:7, groupColumn = 22, order = c(6,2,7,3,5,4,1)) +
  gghighlight::gghighlight(dob == 108 | dob == 95, label_key = dob, use_direct_label = TRUE, keep_scales = TRUE) +
  labs(x = "", y = "") +
  theme +
  theme(legend.position = "none")+
  scale_x_discrete(labels = c("mom", "md", "wd", "mu", "lr", "nl", "bs"))

#bar plot iter 108
#figure 14 left
plot.bar108 = plotShapleyMBO(
  shapley, infill.mbo = "cb", iter.interest = 108, lambda.mbo = 1, 
  sample.size = 20000, type = "bar", avg.phi = FALSE, ci = TRUE, 
  ci.alpha = 0.05, decomp = 1L
)& scale_y_continuous(limits = c(-0.02, 0.003)) # use the same scale as iter 95
# figure 14
plot108 = plot.bar108 + par.plot108


# comparison of lr in iter 108 and 95: all the upcoming stuff put together generates
# table 12 in the thesis

# contributions
contr.lr95 =  c(
  "phi.cb" = shapley95[which(shapley95$feature == "learning_rate"), "phi_cb"] / unique(shapley95$pred.interest_cb - shapley95$pred.average_cb),
  "phi.mean" = shapley95[which(shapley95$feature == "learning_rate"), "phi_mean"] / unique(shapley95$pred.interest_mean - shapley95$pred.average_mean),
  "phi.se" = shapley95[which(shapley95$feature == "learning_rate"), "phi_se"] / unique(shapley95$pred.interest_se - shapley95$pred.average_se)
)
contr.lr108 =  c(
  "phi.cb" = shapley108[which(shapley108$feature == "learning_rate"), "phi_cb"] / unique(shapley108$pred.interest_cb - shapley108$pred.average_cb),
  "phi.mean" = shapley108[which(shapley108$feature == "learning_rate"), "phi_mean"] / unique(shapley108$pred.interest_mean - shapley108$pred.average_mean),
  "phi.se" = shapley108[which(shapley108$feature == "learning_rate"), "phi_se"] / unique(shapley108$pred.interest_se - shapley108$pred.average_se)
)
lr.table = rbind(
  "iter95" = c(
    "lr" = round(opdf[which(opdf$dob == 95), "learning_rate"],4),# lr value original scale
    "avg.lr" = opdf[which(opdf$dob<=108), "learning_rate"] %>% mean() %>% round(digits = 4), # avg lr value upt to that point on original scale
    round(opdf[which(opdf$dob == 95), c("y", "mean", "se", "cb")],4), # error, mean, se and cb predictions
    round(contr.lr95,3) # contributions
  ),
  "iter108" = c(
    "lr" = round(opdf[which(opdf$dob == 108), "learning_rate"],4),
    "avg.lr" = opdf[which(opdf$dob<=108), "learning_rate"] %>% mean() %>% round(digits = 4),
    round(opdf[which(opdf$dob == 108), c("y", "mean", "se", "cb")],4), 
    round(contr.lr108,3)
  )
)
# using opdf.ordered we see that proposals in iter 95 and 108 are respectively the best and second best target values
# achieved in the process
y.order = order(opdf[,"y"], decreasing = FALSE)
opdf.ordered = opdf[y.order,]

# here we compute the empirical distribution function of the mean prediction using the proposals only 
# in order to get the percentile of the mean of proposal 95 and 108
ecdf.m = ecdf(opdf[which(opdf$dob >0), "mean"])
m.95 = ecdf.m(opdf[which(opdf$dob == 95), "mean"])
m.108 = ecdf.m(opdf[which(opdf$dob == 108), "mean"])

# same for the se
ecdf.se = ecdf(opdf[which(opdf$dob >0), "se"])
se.95 = ecdf.se(opdf[which(opdf$dob == 95), "se"])
se.108 = ecdf.se(opdf[which(opdf$dob == 108), "se"])
