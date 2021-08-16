# example for checkSampleSize (see figure 16)
df = data.frame(
  par = rep(c("x1", "x2"),2),
  estimate = rep(c(3,2),2),
  error = c(0, 0.75, 0, 1.25),
  sample.size = rep(c("sufficient K", "increase K"), each = 2)
) %>% tidyr::pivot_longer(cols = c(estimate, error)) %>% dplyr::slice(which(value != 0))



plot = ggplot(df, aes(x = par, y = value, fill = name)) +
  geom_bar(position="stack", stat="identity") +
  labs(y = "phi", x = "", fill = "") +
  facet_grid(cols = vars(sample.size)) +
  scale_fill_manual(values = c("#F8766D", "#999999")) +
  theme_bw()+  
  theme(text = element_text(size = 15))


