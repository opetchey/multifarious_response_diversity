library(ggplot2)

# create data frame
library(truncnorm)

set.seed(46354)
df_low_z1 <- data.frame(E1 = c(rtruncnorm(n=1000, a=250, b=350, mean=270, sd=15), rtruncnorm(n=1000, a=250, b=350, mean=275, sd=15), rtruncnorm(n=1000, a=250, b=350, mean=272, sd=15),
                               rtruncnorm(n=1000, a=50, b=350, mean=280, sd=15)),
                        species = rep(c("Sp 1", "Sp 2", "Sp 3", "Sp 4"), each = 1000))

# plot normal distributions
low_z1_div <- ggplot(df_low_z1, aes(x = E1, fill = species)) +
  geom_density(alpha = 0.5) +
  scale_fill_viridis_d() +
  theme_classic() +
  labs(y = "Growth rate") + ggtitle("Low T optimum diversity") +labs(x = "Temperature", tag = "(g)") +
  theme(
    plot.title = element_text(hjust = 0.5)) +
  theme_classic(base_size = 15) +
  annotate('rect', xmin=250, xmax=305, ymin=0.0, ymax=0.035, alpha=.2, fill='red')+
  annotate('text', x=280, y=0.032, label='Intermediate T mean')+
  annotate('rect', xmin=200, xmax=250, ymin=0.0, ymax=0.035, alpha=.2, fill='green')+
  annotate('text', x=220, y=0.032, label='Low T mean')+
  annotate('rect', xmin=305, xmax=350, ymin=0.0, ymax=0.035, alpha=.2, fill='blue')+
  annotate('text', x=330, y=0.032, label='High T mean')+
  theme(
    plot.title = element_text(hjust = 0.5))

low_z1_div

# create data frame
set.seed(46354)
df_medium_z1 <-data.frame(E1 = c(rtruncnorm(n=1000, a=200, b=285, mean=235, sd=15), rtruncnorm(n=1000, a=210, b=320, mean=265, sd=15), rtruncnorm(n=1000, a=250, b=350, mean=285, sd=15),
                                 rtruncnorm(n=1000, a=250, b=350, mean=288, sd=15)),
                          species = rep(c("Sp 1", "Sp 2", "Sp 3", "Sp 4"), each = 1000))

# plot normal distributions
medium_z1_div <- ggplot(df_medium_z1, aes(x = E1, fill = species)) +
  geom_density(alpha = 0.5) +
  scale_fill_viridis_d() +
  theme_classic() +
  labs(y = "Growth rate") + ggtitle("Medium T optimum diversity") +labs(x = "Temperature", tag = "(h)") +
  theme(
    plot.title = element_text(hjust = 0.5))+
  theme_classic(base_size = 15) +
  annotate('rect', xmin=250, xmax=305, ymin=0.0, ymax=0.035, alpha=.2, fill='red')+
  annotate('text', x=280, y=0.032, label='Intermediate T mean')+
  annotate('rect', xmin=200, xmax=250, ymin=0.0, ymax=0.035, alpha=.2, fill='green')+
  annotate('text', x=220, y=0.032, label='Low T mean')+
  annotate('rect', xmin=305, xmax=350, ymin=0.0, ymax=0.035, alpha=.2, fill='blue')+
  annotate('text', x=330, y=0.032, label='High T mean')+
  theme(
    plot.title = element_text(hjust = 0.5))

medium_z1_div

# create data frame
set.seed(46354)
df_high_z1 <- data.frame(E1 = c(rtruncnorm(n=1000, a=200, b=350, mean=230, sd=15), rtruncnorm(n=1000, a=220, b=300, mean=250, sd=15), rtruncnorm(n=1000, a=250, b=350, mean=285, sd=15),
                                rtruncnorm(n=1000, a=250, b=350, mean=320, sd=15)),
                         species = rep(c("Sp 1", "Sp 2", "Sp 3", "Sp 4"), each = 1000))
# plot normal distributions
high_z1_div <- ggplot(df_high_z1, aes(x = E1, fill = species)) +
  geom_density(alpha = 0.5) +
  scale_fill_viridis_d() +
  theme_classic() +
  labs(y = "Growth rate") + ggtitle("High T optimum diversity") +labs(x = "Temperature", tag = "(i)") +
  theme(
    plot.title = element_text(hjust = 0.5))+
  theme_classic(base_size = 15) +
  annotate('rect', xmin=250, xmax=305, ymin=0.0, ymax=0.035, alpha=.2, fill='red')+
  annotate('text', x=280, y=0.032, label='Intermediate T mean')+
  annotate('rect', xmin=200, xmax=250, ymin=0.0, ymax=0.035, alpha=.2, fill='green')+
  annotate('text', x=220, y=0.032, label='Low T mean')+
  annotate('rect', xmin=305, xmax=350, ymin=0.0, ymax=0.035, alpha=.2, fill='blue')+
  annotate('text', x=330, y=0.032, label='High T mean') +
  theme(
    plot.title = element_text(hjust = 0.5))
high_z1_div

comb <- low_z1_div/ medium_z1_div / high_z1_div & theme(legend.position = "bottom")
fig_optimum <- comb + plot_layout(guides = "collect", heights = c(0.5, 0.5, 0.5))
ggsave(fig_optimum, file="/Users/francesco/Documents/GitHub/multifarious_response_diversity/fig_optimum.jpg", width=8, height=8)

fig_optimum_color <- comb + plot_layout(guides = "collect", heights = c(0.5, 0.5, 0.5))
changing_mean <- (p_E1 + p_E2)/(p2_E1 + p2_E2)/(p3_E1 + p3_E2)
mean_value_env_change <- ggarrange(changing_mean, fig_optimum_color, ncol = 2, nrow = 1)
ggsave(mean_value_env_change, file="/Users/francesco/Documents/GitHub/multifarious_response_diversity/mean_value_env_change.pdf", width=12, height=8)

dd_plot +rdiv_plot + sing_plot 
