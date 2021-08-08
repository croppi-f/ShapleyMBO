library(iml)
library(tidyverse)
library(patchwork)
# common layers for the plots and additional stuff
theme = theme_bw() + 
  theme(text = element_text(size = 15),
        legend.position = "bottom"
  )
hline = geom_hline(yintercept = 0, alpha = 0.5, linetype = "longdash")
###########################################################
###  Discussion Sampling Strategy - Payout Paths       ####
###########################################################
# payout paths for hyper ellipsoid
shapley.he = readRDS(paste0(getwd(),"/analysis/hyper_ellipsoid/shapley1e3.rds"))
payout.he = shapley.he %>% dplyr::select(
  lambda, iter, pred.interest_mean, pred.average_mean, pred.interest_se, pred.average_se, 
  pred.interest_cb, pred.average_cb
) %>% dplyr::slice(seq(1, 15360, 4)) %>% 
  dplyr::group_by(lambda, iter) %>%
  dplyr::summarise(across(.fns = mean))
payout.he.long =  payout.he %>% 
  tidyr::pivot_longer(
    cols = starts_with("pred"),
    names_to = "pred",
    names_ptypes = list(pred = factor()),
    names_prefix = "pred.",
    values_to = "value"
  ) %>%
  tidyr::separate(pred, into = c("pred", "contribution"), sep = "_", remove = TRUE)
# 1) plot actual and the average predictions
plot.payout.he = ggplot(payout.he.long, aes(x = as.numeric(iter), y = value, color = contribution, linetype = pred)) +
  hline +
  geom_line() +
  facet_grid(cols = vars(lambda), rows = vars(contribution), scales = "free_y") +
  scale_color_manual(name = "", values = c("#999999", "#F8766D", "#00BFC4")) +
  scale_linetype_manual(values = c("dashed", "solid"), labels = c("average", "actual")) +
  labs(x = "iter", y = "prediction", linetype = "") +
  theme
#2) plot payout
# scale the payout of se with lambda
payout.he$pred.interest_se_scaled = ifelse(payout.he$lambda == "1", payout.he$pred.interest_se, 10 * payout.he$pred.interest_se)
payout.he$pred.average_se_scaled = ifelse(payout.he$lambda == "1", payout.he$pred.average_se, 10 * payout.he$pred.average_se)
plot.payout.he.ribbon = ggplot(payout.he, aes(x = as.numeric(iter))) +
  hline +
  geom_ribbon(aes(ymin = pred.interest_mean - pred.average_mean, ymax = pred.average_mean - pred.average_mean, fill = "#F8766D"), alpha = 0.5) +
  geom_ribbon(aes(ymin = pred.average_se_scaled - pred.average_se_scaled, ymax = pred.average_se_scaled - pred.interest_se_scaled, fill = "#00BFC4"), alpha = 0.5) +
  scale_fill_identity(name = "", guide = "legend", labels = c("se", "mean")) +
  labs(x = "iter", y = "payout", fill = "") +
  facet_grid(cols = vars(lambda)) +
  theme
plot.payout.he.final = plot.payout.he / plot.payout.he.ribbon +
  plot_annotation(tag_levels = "i")


###########################################################
######  Discussion Phi mean - MLP                  ########
###########################################################
#- here we compare the mean contribution of the proposal in iteration 95 of the 
#  MLP tuning example with the original SM and the final SM to see if results c
#  change
shapley.mlp = readRDS(paste0(getwd(), "/analysis/mlp_phoneme/shapley_phoneme_2e4.rds"))
# take only the results and not the time required
shapley.mlp = shapley.mlp[[1]]
# results with the original sm
shapley.sm95 = shapley.mlp[which(shapley.mlp$iter == "95"), c("feature", "phi_mean")]

# results with the final surrogate model
#- for that we have to explain the istance with iml::Shapley() and it is important
#  that we recreate the same Shapley Object
mbo.mlp = readRDS(paste0(getwd(),"/analysis/mlp_phoneme/mbo_phoneme_lambda_1.rds"))
sm113 = mbo.mlp$models$`113`
opdf.mlp = as.data.frame(mbo.mlp$opt.path)
explicand = opdf.mlp[which(opdf.mlp$dob == 95), c(1:7, 19)]
# reproduce the "same" Shapley Object as ShapleyMBO in iter95
set.seed(1)
samples = ParamHelpers::generateDesign(n = 7000, par.set = mbo.mlp$opt.path$par.set, fun = lhs::randomLHS)
samples$mean = predict(sm113, newdata = samples)$data$response
X = rbind(
  explicand,
  samples, 
  make.row.names = FALSE
)
P = iml::Predictor$new(model = sm113, data = X, y = "mean")
S = iml::Shapley$new(
  predictor = P, x.interest = X[1, ], sample.size = 20000
)
shapley.sm113 = S$results

shapley.all.sm = dplyr::left_join(shapley.sm95, shapley.sm113[, c("feature",  "phi")], by = "feature")
colnames(shapley.all.sm) = c("parameter", "phi_mean_sm95", "phi_mean_sm113")
