library(patchwork)
library(tidyverse)
# common layers for the plots and additional stuff
theme = theme_bw() + 
  theme(text = element_text(size = 15),
        legend.position = "bottom"
  )
lambda.col = c("red", "blue")
hline = geom_hline(yintercept = 0, alpha = 0.5, linetype = "longdash")
par.labels = c(
  expression(theta[1]),
  expression(theta[2]),
  expression(theta[3]),
  expression(theta[4])
)

################################################################################
###############         MBO results           ##################################
################################################################################
opdf_lambda_all = readRDS(paste0(getwd(),"/analysis/hyper_ellipsoid/opdf_lambda_all.rds"))

#1. prediction paths 
#  For each iter and lambda we compute the average
#  over the 30 runs for: y, mean, se and cb
# path is the opdf with dob > 0 
path = dplyr::select(opdf_lambda_all, c(1,6,iter = dob, 11, 16, 17, 18)) %>%
  dplyr::slice(which(iter > 0)) %>% dplyr::group_by(iter, lambda) %>% 
  dplyr::summarise_all(.funs = mean)
# here are the single plots, skip and see directly plot.summary.bo
{
  plot.y = ggplot(path, aes(x = iter, y = y, color = lambda)) +
    geom_line() +
    scale_color_manual(values = lambda.col) +
    theme +
    theme(legend.position = "none")
  plot.y.zoom = plot.y + coord_cartesian(ylim = c(0, 25))
  plot.mean = ggplot(path, aes(x = iter, y = mean, color = lambda)) +
    geom_line() +
    scale_y_continuous(limits = c(min(path$y), max(path$y))) +
    scale_color_manual(values = lambda.col) +
    theme +
    theme(legend.position = "none")
  plot.mean.zoom = plot.mean + coord_cartesian(ylim = c(0, 25))
  plot.se = ggplot(path, aes(x = iter, y = se, color = lambda)) +
    geom_line() +
    scale_color_manual(values = lambda.col) +
    theme +
    theme(legend.position = "none")
  plot.se.zoom = plot.se + coord_cartesian(ylim = c(0, 10))
  plot.cb = ggplot(path, aes(x = iter, y = cb, color = lambda)) +
    geom_line() +
    scale_color_manual(values = lambda.col) +
    theme
  plot.cb.zoom = plot.cb + coord_cartesian(ylim = c(-25, 0)) + 
    theme(legend.position = "none")
  
  }
plot.summary.bo = (plot.y + plot.y.zoom) / 
  (plot.mean + plot.mean.zoom) / 
  (plot.se + plot.se.zoom) /
  (plot.cb + plot.cb.zoom)


# 2. paths of the configurations 
#- dist.to.opt computes the average distance of the proposal to their optimum value,
#  which is 0
# here are the single plots, skip and see directly plot.summary.props
{
  
  dist.to.opt = opdf_lambda_all %>% dplyr::select(iter = dob, lambda, x1,x2,x3,x4) %>%
    dplyr::mutate(x1 = abs(x1), x2 = abs(x2), x3 = abs(x3), x4 = abs(x4)) %>% 
    dplyr::mutate(id = rep(c(1:16, rep(17, 64)), 60)) %>% 
    dplyr::group_by(iter, lambda, id) %>% dplyr::summarise(x1 = mean(x1), x2 = mean(x2), x3 = mean(x3), x4 = mean(x4))
  
  par.plot.props = GGally::ggparcoord(
    dist.to.opt, 
    4:7, 
    groupColumn = 1,
    scale = "globalminmax"
  ) +
    hline +
    scale_color_viridis_c() +
    scale_x_discrete(labels = par.labels) +
    labs(y = "distance to optimum", x = "") +
    facet_grid(cols = vars(lambda)) +
    theme
  
  dist.to.opt2 = opdf_lambda_all %>% dplyr::select(iter = dob, lambda, x1,x2,x3,x4) %>%
    dplyr::slice(which(iter >0)) %>% 
    dplyr::mutate(x1 = abs(x1), x2 = abs(x2), x3 = abs(x3), x4 = abs(x4)) %>% 
    dplyr::group_by(iter, lambda) %>% dplyr::summarise(x1 = mean(x1), x2 = mean(x2), x3 = mean(x3), x4 = mean(x4)) %>% 
    tidyr::pivot_longer(cols = c(x1, x2, x3, x4), names_to = "feature") 
  dist.to.opt2$feature = factor(dist.to.opt2$feature, labels = c("theta[1]", "theta[2]", "theta[3]", "theta[4]"))
  
  plot.props.lw = ggplot(dist.to.opt2, aes(x = iter, y = value, color = feature)) +
    hline +
    geom_line() + 
    scale_color_discrete(labels = par.labels) +
    facet_grid(cols = vars(lambda)) + 
    labs(x = "iter", y = "distance to optimum", color = "") + 
    theme
  plot.props.fw = ggplot(dist.to.opt2, aes(x = iter, y = value, color = lambda)) +
    hline +
    geom_line() + 
    facet_wrap(~feature, labeller = label_parsed) +
    scale_color_manual(values = lambda.col) +
    labs(x = "iter", y = "distance to optimum") + 
    theme
  }
plot.summary.props = plot.props.fw / par.plot.props +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "i") & 
  theme(legend.position = "bottom")
plot.summary.props2 = plot.props.fw / (plot.props.lw / par.plot.props) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "i") & 
  theme(legend.position = "bottom")

################################################################################
###############   Shapley results             ##################################
################################################################################
# load the results
S = readRDS(paste0(getwd(),"/analysis/hyper_ellipsoid/shapley1e3.rds"))
# compute the average contribution and std.dev for each feature, lambda and iter
# average 
avg.shapley = S %>% dplyr::select(iter, feature, lambda, starts_with("phi_")) %>%
  dplyr::group_by(iter, feature,lambda) %>%
  dplyr::summarise_all(.funs = mean)
# pivot to longer format 
avg.shapley.long = avg.shapley %>% dplyr::select(lambda, iter, feature, cb = phi_cb, mean = phi_mean, se = phi_se_scaled) %>%
  tidyr::pivot_longer(cols = c(mean, se, cb), names_to = "contribution", names_ptypes = list(contribution = factor()), values_to = "phi.avg")

# stad.dev
sd.shapley = S %>% dplyr::select(iter, feature, lambda, starts_with("phi_")) %>%
  dplyr::group_by(iter, feature,lambda) %>%
  dplyr::summarise_all(.funs = sd)
sd.shapley.long = sd.shapley %>% dplyr::select(lambda, iter, feature, cb = phi_cb, mean = phi_mean, se = phi_se_scaled) %>%
  tidyr::pivot_longer(cols = c(mean, se, cb), names_to = "contribution", names_ptypes = list(contribution = factor()), values_to = "phi.sd")

#merge the data and use shapley for visualizations 
# feature is transformed in order to display theta symbol in the plots
shapley = dplyr::left_join(avg.shapley, sd.shapley, by = c("iter", "feature", "lambda"), suffix = c(".avg", ".sd"))
shapley$feature = factor(shapley$feature, labels = c("theta[1]", "theta[2]", "theta[3]", "theta[4]"))
# pivot to longer format 
shapley.long = dplyr::left_join(avg.shapley.long, sd.shapley.long, by = c("iter", "feature", "lambda", "contribution"))
shapley.long$feature = factor(shapley.long$feature, labels = c("theta[1]", "theta[2]", "theta[3]", "theta[4]"))

################################################################################
###############   1. Summary table                    ##########################
################################################################################
#- here we provide a summary table of all the contributions over the entire path
contr = shapley.long %>% dplyr::group_by(lambda, feature, contribution) %>% 
  dplyr::summarise(phi = mean(phi.avg), sd = mean(phi.sd))
table.contr = tidyr::pivot_wider(contr, id_cols = c("lambda", "feature"), names_from = "contribution", values_from = c("phi", "sd")) %>%
  dplyr::select(1,2, 5,8, 3,6, 4,7)
table.contr[,3:8] = round(table.contr[,3:8], 2)

#- comparison l10 contribution first 10 and final 10 iterations
contr.l10 = shapley.long %>% dplyr::slice(which(lambda == "10")) %>% 
  dplyr::select(iter, lambda, feature, contribution, phi.avg, phi.sd) %>% 
  dplyr::slice(which((iter %in% as.character(c(1:10, 55:64))))) %>% 
  dplyr::mutate(interval = ifelse(iter %in% as.character(1:10), "begin", "end")) %>% 
  dplyr::group_by(interval, contribution, feature) %>% 
  dplyr::summarise(phi = mean(phi.avg), sd = mean(phi.sd))
table.contr.l10 = tidyr::pivot_wider(contr.l10, id_cols = c("feature", "interval"), names_from = "contribution", values_from = c("phi", "sd")) %>%
  dplyr::select(1,2, 5,8, 3,6, 4,7)
table.contr.l10[,3:8] = round(table.contr.l10[,3:8], 2)  
################################################################################
###############   2. Bar Plots                        ##########################
################################################################################
source("R/plotShapleyHypEllmultiRun.bar.sd.R")
#-here we want the iteration with the most similar proposed configurations between 
# l1 and l10
#-the goal is to assess if there any differences in the contributions between the "same" 
# instance using different lambdas

# prepare the data and find the iter
{
  data.L2.norm = opdf_lambda_all %>% dplyr::select(iter = dob, lambda, x1,x2,x3,x4) %>%
    dplyr::slice(which(iter >0)) %>% 
    dplyr::mutate(x1 = abs(x1), x2 = abs(x2), x3 = abs(x3), x4 = abs(x4)) %>% 
    dplyr::group_by(iter, lambda) %>% dplyr::summarise(x1 = mean(x1), x2 = mean(x2), x3 = mean(x3), x4 = mean(x4)) %>% 
    dplyr::select(iter, x1, x2, x3, x4)
  # compute the L2 norm between points
  L2.norm = list()
  for(i in 1:64) {
    data = data.L2.norm[which(data.L2.norm$iter == i), ]
    L2.norm[[i]] = dist(data)
    names(L2.norm)[i] = i
  }
  L2.norm = unlist(L2.norm) 
  min.L2 = which.min(L2.norm) # 59
}
#- im iter 59 instances are the most similar
bar59.l1 =  plotShapleyHypEllmultiRun.bar.sd(S, lambda.mbo = 1, iter.interest = 59)
bar59.l10 =  plotShapleyHypEllmultiRun.bar.sd(S, lambda.mbo = 10, iter.interest = 59)
#- we select only a subset of plots. Notice that since patchwork is used single
#  plots can be accesses like list elements
bar59 = (bar59.l1[[1]] + bar59.l1[[2]]) / 
  (bar59.l10[[1]] + bar59.l10[[2]]) +
  plot_layout(guides = "collect") &
  plot_annotation(tag_levels = "i") & 
  theme(text = element_text(size = 15), legend.position = "bottom") & 
  scale_y_continuous(limits = c(-45, 21))

#summary tables for iter 59
#prepare the data (named shapley59)
{
  # average contribution
  shapley59 = S %>% dplyr::slice(which(iter == "59")) %>% dplyr::select(
    lambda, feature, phi_mean, pred.interest_mean, pred.average_mean,
    phi_se_scaled, pred.interest_se, pred.average_se,
    phi_cb, pred.interest_cb, pred.average_cb
  )  %>% group_by(lambda, feature) %>%
    dplyr::summarise(across(c(starts_with("phi"), starts_with("pred")), mean))
  shapley59$weight.phi_cb = round(abs(shapley59$phi_cb / (shapley59$pred.interest_cb - shapley59$pred.average_cb)), 4)
  shapley59$weight.phi_mean = round(abs(shapley59$phi_mean / (shapley59$pred.interest_mean - shapley59$pred.average_mean)), 4)
  shapley59$weight.phi_se = ifelse(
    shapley59$lambda == "1",
    round(abs(shapley59$phi_se_scaled / (shapley59$pred.interest_se - shapley59$pred.average_se)), 4),
    round(abs(shapley59$phi_se_scaled / (10 * shapley59$pred.interest_se - 10 * shapley59$pred.average_se)), 4)
  )
  # stad.dev
  shapley59.sd = S %>% dplyr::slice(which(iter == "59")) %>% dplyr::select(
    lambda, feature, phi_mean.sd = phi_mean, phi_se_scaled.sd = phi_se_scaled, phi_cb.sd = phi_cb
  )  %>% group_by(lambda, feature) %>%
    dplyr::summarise(across(starts_with("phi"), sd))
  # merge the data
  shapley59 = dplyr::left_join(shapley59, shapley59.sd, by = c("lambda", "feature"))
}
table59.pred = round(
  unique(shapley59[, c(10,11,6,7,8,9)]),
  2
)
#for cb, mean, se: phi, sd, weight.phi
table59.phi = round(
  rbind(
    l1 = t(shapley59[1:4, c(5,17,12, 3,15,13, 4,16,14)]),
    l10 = t(shapley59[5:8, c(5,17,12, 3,15,13, 4,16,14)])
  ),
  2
)
  
# exploration factor: how much have the the dims explored until iter 58
#- prepare data: with and without init design we compute the std.dev of the config
# values as a measure of exploretaion
{
  explored.init.des = dist.to.opt %>% dplyr::slice(which(id <= 16)) %>% dplyr::select(-id)%>% 
    tidyr::pivot_longer(cols = c(x1, x2, x3, x4), names_to = "feature") 
  explored.init.des $feature = factor(explored.init.des $feature, labels = c("theta[1]", "theta[2]", "theta[3]", "theta[4]"))
  explored.with.id = dist.to.opt2 %>% dplyr::slice(which(iter <= 58)) %>%
    dplyr::bind_rows(explored.init.des) %>% 
    dplyr::group_by(feature, lambda) %>% 
    dplyr::summarise(sd = sd(value)) %>% tidyr::pivot_wider(names_from = feature, values_from = sd)
  explored.no.id = dist.to.opt2 %>% dplyr::slice(which(iter <= 58)) %>%
    dplyr::group_by(feature, lambda) %>% 
    dplyr::summarise(sd = sd(value)) %>% tidyr::pivot_wider(names_from = feature, values_from = sd) 
}
explored = dplyr::bind_rows(explored.with.id, explored.no.id, .id = "with.id")%>%
  dplyr::mutate(with.id = rep(c("yes", "no"),1, each = 2))
explored[,3:6] = round(explored[,3:6], 2)


################################################################################
###############   3. Desirability paths               ##########################
################################################################################
#-mid.<bla> is used to measure the average contributions (over the entire path) of the pars
# in the lambda1 scenario and display it in the lambda10 plot for a comparison
{
  mid.cb.l1 = geom_hline(
    data = cbind(
      contr[which(contr$contribution == "cb" & contr$lambda == "1"), c(2,4)],
      lambda = "10"
    ),
    mapping = aes(group = feature, yintercept = phi),
    linetype = "dotdash"
  )
  mid.m.l1 = geom_hline(
    data = cbind(
      contr[which(contr$contribution == "mean" & contr$lambda == "1"), c(2,4)],
      lambda = "10"
    ),
    mapping = aes(group = feature, yintercept = phi),
    linetype = "dotdash"
  )
  mid.se.l1 = geom_hline(
    data = cbind(
      contr[which(contr$contribution == "se" & contr$lambda == "1"), c(2,4)],
      lambda = "10"
    ),
    mapping = aes(group = feature, yintercept = phi),
    linetype = "dotdash"
  )
}
# overall desirability: phi(cb)
phi.cb.path = ggplot(data = shapley, aes(x = as.numeric(iter), y = phi_cb.avg)) +
  hline +
  geom_line(size = 0.75, color = "#999999") +
  geom_errorbar(aes(ymin = phi_cb.avg - phi_cb.sd, ymax = phi_cb.avg + phi_cb.sd), color = "#999999") +
  mid.cb.l1 +
  scale_color_discrete(labels = par.labels) +
  labs(x = "iter", y = "phi (cb)") +
  facet_grid(cols = vars(feature), rows = vars(lambda), labeller = label_parsed) +
  theme

# explore-exploit trade-off (phi_mean and phi_se paths)
eeto.path = ggplot(data = shapley.long[which(shapley.long$contribution != "cb"), ], aes(x = as.numeric(iter), y = phi.avg, color = contribution)) +
  hline +
  geom_line(size = 0.75) +
  geom_errorbar(aes(ymin = phi.avg - phi.sd, ymax = phi.avg + phi.sd)) +
  mid.m.l1 +
  mid.se.l1 +
  labs(x = "iter", y = "phi", color = "contribution") +
  facet_grid(cols = vars(feature), rows = vars(lambda), labeller = label_parsed) +
  theme

#-here we investigate further the results of the mean and se paths
#-we take into consideration the first 10 iters of the process
#1. display the same parallel plot as above but highlight the first 10 iters 
#  (notice decreasing trend --> lower pars more explored in initial iters and this
#   why their se contr increase in the beginning)
par.plot.props.10dob = GGally::ggparcoord(
  dist.to.opt[which(dist.to.opt$lambda == "10"), ], 
  4:7, 
  groupColumn = 1,
  scale = "globalminmax"
) +
  scale_color_viridis_c() + 
  scale_x_discrete(labels = par.labels)+
  gghighlight::gghighlight(iter > 0 & iter <= 5, use_direct_label = FALSE, keep_scales = TRUE) +
  hline +
  labs(y = "distance to optimum", x = "") +
  theme

#2. histogram for the high variance in the contributions of lambda10
hist.l10 = opdf_lambda_all[which(opdf_lambda_all$dob > 0 & opdf_lambda_all$dob <= 10 & opdf_lambda_all$lambda == "10"), c(2:5, 7)] %>% 
  tidyr::pivot_longer(1:4, names_to = "feature")
hist.l10$feature = factor(hist.l10$feature, labels = c("theta[1]", "theta[2]", "theta[3]", "theta[4]"))

plot.hist.l10 = ggplot(hist.l10, aes(x = value), fill = "#999999") +
  geom_histogram(bins = 15, color = "black") +
  facet_grid(cols = vars(feature), labeller = label_parsed) +
  labs(x = "value", y = "count", fill = "") + 
  theme
data.hist = data.frame(
  parameter = rep(c("theta[1]", "theta[2]", "theta[3]", "theta[4]"), each = 15),
  count = as.numeric(ggplot_build(plot.hist.l10)$data[[1]]$count)
)


