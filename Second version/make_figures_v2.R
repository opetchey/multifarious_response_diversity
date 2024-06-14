
### Figures Second version
# Manuscript figure

# rm(list = ls())

# packages
pkgs <- c("here", "dplyr", "ggplot2", "mgcv", "gratia", "scales", "tidyverse", "purrr",
          "gridExtra", "cowplot", "ggsci", "viridis", "pracma")
vapply(pkgs, library, logical(1L), logical.return = TRUE, character.only = TRUE)

source.files <- list.files(here("r"), full.names = TRUE)
sapply(source.files, source, .GlobalEnv)

load(here("Data/my_work_space.RData"))
# the above loads all these old gratia functions - you should use the dev
# version of gratia and install the binary from my r-universe
# Now I delete all the old functions that are stored in the sourced image
rm(list = ls(getNamespace("gratia"), all.names = TRUE))

############ Figure 1 ################
# Fig 1



E1_series <- seq(0, 10, 0.2)
E2_series <- seq(0, 10, 0.2)
# Simulate spp performance curves with the modified Eppley function with and without interactive effect.

# Without interaction
set.seed(2465)
s <- 4 ## number of species
a1_mean <- 1e-9
b1_mean <- 0.063
z1_mean <- 5
w1_mean <- 3
a2_mean <- 1e-3
b2_mean <- 0.02
z2_mean <- 5
w2_mean <- 2
zint_mean <- 0 ## no interaction
sd_rate_mean <- 0.02
z1_range <- 5
z2_range <- 3
par_table <- tibble(a1 = rnorm(s, a1_mean, 0),#a1_mean/100),
                    b1 = rnorm(s, b1_mean, 0),#b1_mean/100),
                    z1 = runif(s, min = z1_mean - z1_range/2, max = z1_mean + z1_range / 2),
                    w1 = rnorm(s, w1_mean, 0),#w1_mean/100),
                    a2 = rnorm(s, a2_mean, 0),#a2_mean/100),
                    b2 = rnorm(s, b2_mean, 0),#b2_mean/100),
                    z2 = runif(s, min = z2_mean - z2_range/2, max = z2_mean + z2_range / 2),
                    w2 = rnorm(s, w2_mean, 0),#w2_mean/100),
                    z_int21 = rnorm(s, zint_mean, zint_mean/10),
                    sd_rate = rnorm(s, sd_rate_mean, 0)
)


## The normal house keeping
par_list <- Partable_2_parlist(par_table)
species_pars_no_inter <- tibble(species = paste0("s", 1:s), pars = par_list) %>%
  rowwise() %>%
  mutate(expt = Make_expt(E1_series, E2_series, pars))



species_pars_no_inter <- species_pars_no_inter %>% unnest(expt)
species_pars_no_inter %>%  group_by(species) %>%  nest() 


(dd1 <- species_pars_no_inter %>% dplyr::select(species, E1, E2, rate))

#Fit response surface for each sp (done with GAMs)

# Create new env data
new_data <- tibble(E1 = seq(0, 10, 0.1),
                   E2 = seq(0, 10, 0.1))

# Without interaction
# GAMs
nested_gams_no_inter <- species_pars_no_inter %>% 
  nest(cols =-species) %>% 
  mutate(
    gams = map(cols, ~ gam(rate ~ te(E1, E2, k = c(10, 10)),
                           data = .x,
                           method = "REML")),
    predicted = map(gams, ~ predict(.x, newdata = new_data))
  )



# Get gams coefficients
nested_gams_no_inter %>% 
  mutate(coefs = map(gams, tidy, conf.int = TRUE)) %>% 
  unnest(coefs)


# Get gams glance
nested_gams_no_inter %>% 
  mutate(results = map(gams, glance), 
         R.square = map_dbl(gams, ~ summary(.)$r.sq))  %>% 
  unnest(results) 


sp1 <- draw(nested_gams_no_inter$gams[[1]], rug = FALSE) &
  labs(tag = "(a)",  caption = "")
sp2 <- draw(nested_gams_no_inter$gams[[2]], rug = FALSE) &
  labs(tag = "(b)",  caption = "")
sp3 <- draw(nested_gams_no_inter$gams[[3]], rug = FALSE) &
  labs(tag = "(c)",  caption = "")
sp4 <- draw(nested_gams_no_inter$gams[[4]], rug = FALSE) &
  labs(tag = "(d)",  caption = "")

#Data wrangling and partials derivatives calculations

# creating new environmental conditions 
refs2 <- crossing(E1 = seq(0, 10, 0.2),
                  E2 = seq(0, 10, 0.2))

#community 1 
m_list_no <- (nested_gams_no_inter$gams)

my_spp_names_no <- (nested_gams_no_inter$species)


(pd_list_no <- modify_depth(m_list_no, 1, ~ get_partials(., refs2)))




# create a dataframe with E1 and E2 changing over time
set.seed(64)

library(MASS)


# Generate multivariate normal data without constraints
refs_ts <- mvrnorm(n = 10, 
                   mu = c(mean(E1_series), mean(E2_series)),
                   Sigma = matrix(c(9, 0.0, 0.0, 9), nrow = 2),
                   empirical = TRUE,  # Ensure data is distributed as specified
                   tol = 1e-6)         # Tolerance for convergence



# refs_ts <- mvrnorm(n = 10,        #  specifies the sample size
#                    mu = c(mean(E1_series), mean(E2_series)), #specifies the mean values of each column
#                    Sigma = matrix(c( 9, 0.0,
#                                      0.0, 9), #  specifies the correlation matrix
#                                   nrow = 2))
# 

refs_ts <- refs_ts %>%  as_tibble(refs_ts) %>% 
  dplyr::rename(E1 = "V1", E2 = "V2")

refs_ts <- refs_ts%>% 
  dplyr::mutate(time = seq.int(nrow(refs_ts)))

p_E1 <- refs_ts %>% 
  ggplot(aes(x = time, y = E1)) + geom_line(color = "orange", size = 1.5) +
  geom_point(size = 3)+
  theme_classic(base_size = 20) +
  labs(y = "E1", x = "Time", tag = "(a)") + 
  scale_x_continuous(breaks=seq(0, 10, 1))+
  scale_y_continuous(breaks=seq(0, 20, 2))

p_E2 <- refs_ts %>% 
  ggplot(aes(x = time, y = E2)) + geom_line(color = "blue", size = 1.5) +
  geom_point(size = 3)+
  theme_classic(base_size = 20) +
  labs(y = "E2", x = "Time",tag = "(b)") + 
  scale_x_continuous(breaks=seq(0, 10, 1)) +
  scale_y_continuous(breaks=seq(0, 20, 2))

# Plot environmental change over time
p_E1 + p_E2
indipendent_fluctuations <- p_E1 + p_E2

detach("package:MASS", unload=TRUE)


### Response surfaces with change in environmental conditions

sp1_nocor <- sp1 + geom_point(data = refs_ts, mapping = aes(x = E1, y = E2), size = 3) +
  geom_path(data = refs_ts, aes(color = "#0072B2"), size = 1.5) +
  geom_label(data = refs_ts, mapping = aes(label = time)) +
  guides(color="none" , fill=guide_legend(title="Growth \n Rate"))

sp1_nocor <- sp1 + geom_point(data = refs_ts, mapping = aes(x = E1, y = E2), size = 3) +
  geom_label(data = refs_ts, mapping = aes(label = time)) +
  geom_path(data = refs_ts, aes(color = time), arrow = arrow())+ guides(color="none")

sp2_nocor <- sp2 + geom_point(data = refs_ts, mapping = aes(x = E1, y = E2), size = 3) +
  geom_label(data = refs_ts, mapping = aes(label = time)) +
  geom_path(data = refs_ts, aes(color = time), arrow = arrow())+ guides(color="none")

sp3_nocor <- sp3 + geom_point(data = refs_ts, mapping = aes(x = E1, y = E2), size = 3) +
  geom_label(data = refs_ts, mapping = aes(label = time)) +
  geom_path(data = refs_ts, aes(color = time), arrow = arrow())+ guides(color="none")

sp4_nocor <- sp4 + geom_point(data = refs_ts, mapping = aes(x = E1, y = E2), size = 3) +
  geom_label(data = refs_ts, mapping = aes(label = time)) +
  geom_path(data = refs_ts, aes(color = time), arrow = arrow())+ guides(color="none")

ggarrange(sp1_nocor, sp2_nocor, sp3_nocor,  sp4_nocor, ncol=4, nrow=1, common.legend = TRUE, legend="right")

# get partial derivatives - measured values
(pd_list_no <- modify_depth(m_list_no, 1, ~ get_partials(., refs_ts)))
# from list to tibble
(pd_spp_no <- tibble(
  E1_ref = map(pd_list_no, "E1"),
  E2_ref = map(pd_list_no, "E2"),
  pd_E1 = map(pd_list_no, "pd_E1"),
  pd_E2 = map(pd_list_no, "pd_E2")) %>%
    dplyr::mutate(sp = my_spp_names_no) %>%
    relocate(sp, E1_ref, E2_ref, pd_E1, pd_E2) %>%
    unnest(E1_ref, E2_ref, pd_E1, pd_E2))
# add time
(pd_spp_no <-cbind(pd_spp_no, refs_ts$time) %>% 
    dplyr::rename(time = "refs_ts$time"))

# calculation next value for directional derivatives, and get directional derivatives
(pd_spp_no <- pd_spp_no %>% transform( nxt_value_E1 = c(E1_ref[-1], NA)) %>%
    transform(nxt_value_E2 = c(E2_ref[-1], NA)) %>%
    dplyr::mutate(del_E1 = nxt_value_E1 - E1_ref,
                  del_E2 = nxt_value_E2 - E2_ref,
                  unit_vec_mag =  sqrt(del_E1^2 + del_E2^2),
                  uv_E1 = del_E1 / unit_vec_mag,
                  uv_E2 = del_E2 / unit_vec_mag,
                  dir_deriv = pd_E1 * uv_E1 +  pd_E2 * uv_E2) %>% 
    filter(time != max(time)))



# community 1 
# reduce the dataframe and keep only what we need
red_spp <- pd_spp_no %>%  dplyr::select(sp, time, E1_ref, E2_ref, dir_deriv)

red_spp_s1 <- red_spp %>% dplyr::filter(sp == "s1")

dir_plot <- ggplot()+  geom_line(data = red_spp_s1, aes(x = time + 0.5, y = dir_deriv), color = "black", size = 1) +
  geom_point(data = red_spp_s1, aes(x = time + 0.5, y = dir_deriv), color = "black", size = 3) +
  scale_x_continuous(breaks = seq(0, 10, 1), expand = c(0, 0.5)) +
  geom_point(size = 3)+
  labs(x = "Time",y = "Directional derivative", tag = "(d)") +
  geom_hline(yintercept = 0, linetype= "dashed") + theme_classic(base_size = 15) +
  scale_x_continuous(breaks=seq(0, 10, 1))


# View combined plot
dir_plot

# Sp response surface

sp_plot <- sp1 + geom_point(data = refs_ts, mapping = aes(x = E1, y = E2), size = 3) +
  geom_path(data = refs_ts, size = 1.5) +
  geom_point(data = refs_ts, mapping = aes(x = E1, y = E2), size = 1) +
  geom_path(data = refs_ts,  size = 0.5, alpha = 0.5) + guides(color="none" , fill=guide_colourbar(title="Growth \n Rate")) +
  geom_label(data = refs_ts, mapping = aes(label = time), size = 10) + labs(title = "", tag = "(c)")+ theme_classic(base_size = 20) #+ ylab("Salinity (ppt)") + xlab("Temperature (k)")  + 



Fig1 <- (p_E1 + p_E2 + sp_plot) / dir_plot

ggsave(Fig1, file = here("Second version/figures_v2/Fig1.jpg"), width = 22, height = 15)
ggsave(Fig1, file = here("Second version/figures_v2/Fig1.png"), width = 22, height = 15)
ggsave(Fig1, file = here("Second version/figures_v2/Fig1.pdf"), width = 22, height = 15)


###### Figure 2 ############
# This creates the base, then was edited in power point
get_partials <- function(m, refs) {
  refs$pd_x <- NA
  refs$pd_z <- NA
  for(i in 1:nrow(refs)) {
    refs$pd_x[i] <- partial_derivatives(m,
                                        data = refs[i,],
                                        type = "central",
                                        focal = "x")$partial_deriv
    refs$pd_z[i] <- partial_derivatives(m,
                                        data = refs[i,],
                                        type = "central",
                                        focal = "z")$partial_deriv
  }
  refs
}


## Species 1 

df <- data_sim("eg2", n = 2000, dist = "normal", scale = 0.5, seed = 55)

# fit the GAM (note: for execution time reasons, k is set articifially low)
m <- gam(y ~ te(x, z, k = c(5, 5)), data = df, method = "REML")

# draw te(x,z)
p1 <- draw(m, rug = FALSE,
           n_contour = 30,
           contour_col = "lightgrey") & 
  labs(title = "", fill = "Growth", caption = "")

refs <- crossing(x = seq(0.1, 0.9, 0.1),
                 z = seq(0.1, 0.9, 0.1))

refs <- get_partials(m, refs) %>% 
  dplyr::mutate(x_ref = 0.5, z_ref = 0.5) 

radius <- 0.15
num_arrows <- 20
dd1 <- tibble(angle = rep(seq(0, 2*pi, length = num_arrows), nrow(refs))) %>%
  mutate(x_ref = rep(refs$x_ref, each = num_arrows),
         z_ref = rep(refs$z_ref, each = num_arrows)) %>%
  full_join(refs) %>%
  mutate(x = cos(angle) * radius,
         z = sin(angle) * radius,
         dir_deriv = x * pd_x + z * pd_z,
         unit_vec_mag = sqrt(x^2 + z^2))
fi2_a <- p1 +
  geom_segment(data = dd1,
               aes(x = x_ref, y = z_ref,
                   xend = x_ref + x, yend = z_ref + z,
                   col = dir_deriv), alpha = 0.7,
               size=1) +
  scale_colour_gradient2(low = muted("blue"),
                         mid = "white",
                         high = muted("red"),
                         midpoint = 0)


## Fig 2 a - one point

fig2_a <- fi2_a+ labs(tag = "(a)") + guides(fill=guide_colourbar(title="Growth \n Rate")) +
  theme(legend.key.width = unit(0.28, "cm"),  # Adjust width of the legend key
        legend.key.height = unit(0.28, "cm")) + 
  labs(color = "Directional \n  Derivative", title = "Sp 1", , x = "E1", y = "E2") +
  theme_classic(base_size = 20) 

## Fig 2 c - grid of points


p1.2 <- draw(m, rug = FALSE,
             n_contour = 10,
             contour_col = "black")&
  scale_fill_gradient(low = "white", 
                      high = "darkgrey") & 
  labs(title = "", fill = "Growth", caption = "") & theme_classic(base_size = 20) 

refs <- crossing(x = seq(0.1, 0.9, 0.1),
                 z = seq(0.1, 0.9, 0.1))

refs <- get_partials(m, refs) %>% 
  dplyr::rename(x_ref = x, z_ref = z) 

radius <- 0.05
num_arrows <- 20
dd1 <- tibble(angle = rep(seq(0, 2*pi, length = num_arrows), nrow(refs))) %>%
  mutate(x_ref = rep(refs$x_ref, each = num_arrows),
         z_ref = rep(refs$z_ref, each = num_arrows)) %>%
  full_join(refs) %>%
  mutate(x = cos(angle) * radius,
         z = sin(angle) * radius,
         dir_deriv = x * pd_x + z * pd_z,
         unit_vec_mag = sqrt(x^2 + z^2))

fig2_c <- p1.2 +
  geom_segment(data = dd1,
               aes(x = x_ref, y = z_ref,
                   xend = x_ref + x, yend = z_ref + z,
                   col = dir_deriv), alpha = 0.7,
               size=1) +
  scale_colour_gradient2(low = muted("blue"),
                         mid = "white",
                         high = muted("red"),
                         midpoint = 0) + 
  labs(x = "E1", y = "E2", tag = "(c)", title = "Sp 1")+
  theme_classic(base_size = 20) 


### Fig 2e
refs <- crossing(x = seq(0.1, 0.9, 0.35),
                 z = seq(0.1, 0.9, 0.35))
refs <- get_partials(m, refs) %>%
  dplyr::rename(x_ref = x, z_ref = z)
radius <- 0.15
num_arrows <- 4
dd1 <- tibble(angle = rep(seq(0, 2*pi, length = num_arrows), nrow(refs))) %>%
  mutate(x_ref = rep(refs$x_ref, each = num_arrows),
         z_ref = rep(refs$z_ref, each = num_arrows)) %>%
  full_join(refs) %>%
  mutate(x = cos(angle) * radius,
         z = sin(angle) * radius,
         dir_deriv = x * pd_x + z * pd_x,
         unit_vec_mag = sqrt(x^2 + z^2))
fig2_e<- p1.2 +
  geom_segment(data = dd1,
               aes(x = x_ref, y = z_ref,
                   xend = x_ref + x, yend = z_ref + z,
                   col = dir_deriv), alpha = 0.7,
               size=3) +
  scale_colour_gradient2(low = muted("blue"),
                         mid = "white",
                         high = muted("red"),
                         midpoint = 0) +
  geom_point(data = dd1,
             aes(x = x_ref, y = z_ref,
                 size = 2))+
  xlab("E1") + ylab("E2") + 
  labs(title = "Sp 1", tag = "(e)", caption = "") + theme_classic(base_size = 20) +
  guides(size = "none") 


#### Species 2 
### Fig 2 b 
df2 <- data_sim("eg2", n = 2000, dist = "normal", scale = 0.9, seed = 55)

# fit the GAM (note: for execution time reasons, k is set articifially low)
m2 <- gam(y ~ te(x, z, k = c(5, 5)), data = df2, method = "REML")

# draw te(x,z)
p2 <- draw(m2, rug = FALSE,
           n_contour = 30,
           contour_col = "lightgrey") & 
  labs(title = "", fill = "Growth", caption = "")

ref <- crossing(x = seq(0.1, 0.9, 0.1),
                z = seq(0.1, 0.9, 0.1))

ref <- get_partials(m2, ref) %>% 
  dplyr::mutate(x_ref = 0.5, z_ref = 0.5) 

radius <- 0.15
num_arrows <- 20
dd1 <- tibble(angle = rep(seq(0, 2*pi, length = num_arrows), nrow(ref))) %>%
  mutate(x_ref = rep(ref$x_ref, each = num_arrows),
         z_ref = rep(ref$z_ref, each = num_arrows)) %>%
  full_join(ref) %>%
  mutate(x = cos(angle) * radius,
         z = sin(angle) * radius,
         dir_deriv = x * pd_x + z * pd_z,
         unit_vec_mag = sqrt(x^2 + z^2))
fig2_b <- p2 +
  geom_segment(data = dd1,
               aes(x = x_ref, y = z_ref,
                   xend = x_ref + x, yend = z_ref + z,
                   col = dir_deriv), alpha = 0.7,
               size=1) +
  scale_colour_gradient2(low = muted("blue"),
                         mid = "white",
                         high = muted("red"),
                         midpoint = 0)


## Fig 2 b - one point

fig2_b <- fig2_b+ labs(tag = "(b)") + guides(fill=guide_colourbar(title="Growth \n Rate")) +
  theme(legend.key.width = unit(0.28, "cm"),  # Adjust width of the legend key
        legend.key.height = unit(0.28, "cm")) + 
  labs(color = "Directional \n  Derivative", title = "Sp 2", x = "E1", y = "E2") +
  theme_classic(base_size = 20) 

## Fig 2 d - grid of points


p2.2 <- draw(m2, rug = FALSE,
             n_contour = 10,
             contour_col = "black")&
  scale_fill_gradient(low = "white", 
                      high = "darkgrey") & 
  labs(title = "", fill = "Growth", caption = "") & theme_classic(base_size = 20) 

refs <- crossing(x = seq(0.1, 0.9, 0.1),
                 z = seq(0.1, 0.9, 0.1))

refs <- get_partials(m2, refs) %>% 
  dplyr::rename(x_ref = x, z_ref = z) 

radius <- 0.05
num_arrows <- 20
dd1 <- tibble(angle = rep(seq(0, 2*pi, length = num_arrows), nrow(refs))) %>%
  mutate(x_ref = rep(refs$x_ref, each = num_arrows),
         z_ref = rep(refs$z_ref, each = num_arrows)) %>%
  full_join(refs) %>%
  mutate(x = cos(angle) * radius,
         z = sin(angle) * radius,
         dir_deriv = x * pd_x + z * pd_z,
         unit_vec_mag = sqrt(x^2 + z^2))

fig2_d <- p2.2 +
  geom_segment(data = dd1,
               aes(x = x_ref, y = z_ref,
                   xend = x_ref + x, yend = z_ref + z,
                   col = dir_deriv), alpha = 0.7,
               size=1) +
  scale_colour_gradient2(low = muted("blue"),
                         mid = "white",
                         high = muted("red"),
                         midpoint = 0) + 
  labs(x = "E1", y = "E2", tag = "(d)", title = "Sp 2")+
  theme_classic(base_size = 20) 


### Fig 2f
refs <- crossing(x = seq(0.1, 0.9, 0.35),
                 z = seq(0.1, 0.9, 0.35))
refs <- get_partials(m2, refs) %>%
  dplyr::rename(x_ref = x, z_ref = z)
radius <- 0.15
num_arrows <- 4
dd1 <- tibble(angle = rep(seq(0, 2*pi, length = num_arrows), nrow(refs))) %>%
  mutate(x_ref = rep(refs$x_ref, each = num_arrows),
         z_ref = rep(refs$z_ref, each = num_arrows)) %>%
  full_join(refs) %>%
  mutate(x = cos(angle) * radius,
         z = sin(angle) * radius,
         dir_deriv = x * pd_x + z * pd_x,
         unit_vec_mag = sqrt(x^2 + z^2))
fig2_f<- p1.2 +
  geom_segment(data = dd1,
               aes(x = x_ref, y = z_ref,
                   xend = x_ref + x, yend = z_ref + z,
                   col = dir_deriv), alpha = 0.7,
               size=3) +
  scale_colour_gradient2(low = muted("blue"),
                         mid = "white",
                         high = muted("red"),
                         midpoint = 0) +
  geom_point(data = dd1,
             aes(x = x_ref, y = z_ref,
                 size = 2))+
  xlab("E1") + ylab("E2") + 
  labs(title = "Sp 2", tag = "(f)", caption = "") + theme_classic(base_size = 20) +
  guides(size = "none") 




# Combine the plots into a single patchwork layout with three rows and two columns
combined_plots <- wrap_plots(
  fig2_a + theme(legend.position = "none"),
  fig2_b + theme(legend.position = "none"),
  fig2_c + theme(legend.position = "none"),
  fig2_d + theme(legend.position = "none"),
  fig2_e + theme(legend.position = "none"),
  fig2_f + theme(legend.position = "none"),
  nrow = 3
)

legend <- get_legend(
  fig2_b + theme(
    legend.title = element_text(size = 15),  # Increase legend title size
    legend.text = element_text(size = 12),   # Increase legend text size
    legend.key.size = unit(1.5, "lines")     # Increase legend key size
  )
)

# Combine the plots and the legend in a grid layout, adjusting spacing
final_plot <- plot_grid(
  combined_plots,
  legend,
  ncol = 2,
  rel_widths = c(3, 0.5),  # Adjust space between plots and legend
  align = 'v'
) +
  theme(plot.margin = margin(0, 0, 0, 0))  # Adjust margins to reduce white space

# Display the final plot
print(final_plot)


ggsave(final_plot, file = here("Second version/figures_v2/Fig2.png"), width = 22, height = 22)

ggsave(fig2a, file = here("Second version/figures_v2/Fig2a.png"), width = 22, height = 15)
ggsave(fig2b, file = here("Second version/figures_v2/Fig2b.png"), width = 22, height = 15)
ggsave(fig2d, file = here("Second version/figures_v2/Fig2d.png"), width = 22, height = 15)
ggsave(fig2e, file = here("Second version/figures_v2/Fig2e.png"), width = 22, height = 15)



########## Figure 3 #################
# Plotting Response Capacity over the surface 
diverging_palette <- viridis(100, option = "plasma")
### Community 1
#### Create community
E1_series <- seq(0, 50, 1)
E2_series <- seq(0, 50, 1)

## define new_data for following GAM fitting

new_data <- tibble(E1 =  seq(0, 50, 0.2),
                   E2 = seq(0, 50, 0.2))
## Set species parameters
species_pars <- list(a1_mean = 1e-3,
                     b1_mean = 0.02,
                     z1_mean = 20,
                     w1_mean = 10,
                     a2_mean = 1e-3,
                     b2_mean = 0.02,
                     z2_mean = 20,
                     w2_mean = 10,
                     zint_mean = 0,
                     sd_rate_mean = 0,
                     z1_range = NA,
                     z2_range = NA)


species_pars$z1_range <- 50
species_pars$z2_range <- 35

# create reproducible community
s <- 10
set.seed(2456)
comm <- get_community(s, species_pars)

# Visualise spp performances 
performance_comm3 <- multi_performance_plot(comm)
performance_comm3


refs <- crossing(E1 = seq(0, 50, 5),
                 E2 = seq(0, 50, 5))

nested_gams <- comm %>% 
  nest(cols =-species) %>% 
  mutate(
    gams = map(cols, ~ gam(rate ~  te(E1, E2),
                           data = .x,
                           method = "REML")),
    predicted = map(gams, ~ predict(.x, newdata = new_data))
  )

#list of gams
m_list <- (nested_gams$gams)
#list of spp names
my_spp_names <- (nested_gams$species)
# get partial derivatives
pd_list <- modify_depth(m_list, 1, ~ get_partials(., refs))

# from list to tibble
pd_spp <- tibble(
  E1_ref = map(pd_list, "E1"),
  E2_ref = map(pd_list, "E2"),
  pd_E1 = map(pd_list, "pd_E1"),
  pd_E2 = map(pd_list, "pd_E2")) %>%
  dplyr::mutate(sp = my_spp_names) %>%
  relocate(sp, E1_ref, E2_ref, pd_E1, pd_E2) %>%
  unnest(E1_ref, E2_ref, pd_E1, pd_E2)


radius <- 1
num_arrows <- 1000
pd_spp <- crossing(angle = rep(seq(0, 2*pi, length = num_arrows)),
                   E1_ref = pd_spp$E1_ref,
                   E2_ref =pd_spp$E2_ref) %>%
  full_join(pd_spp) %>%
  mutate(E1 = cos(angle) * radius,
         E2 = sin(angle) * radius,
         dir_deriv = E1 * pd_E1 + E2 * pd_E2,
         unit_vec_mag = sqrt(E1^2 + E2^2))


# filter useless stuff
pd_spp <- pd_spp %>% 
  dplyr::select(angle, sp, E1_ref, E2_ref, dir_deriv)

#
# Calculating diversity for each direction and location this should be always the first step
# Dissimilarity
Div_loc_dir <- pd_spp %>% 
  dplyr::group_by(E1_ref, E2_ref, angle) %>% 
  summarise(div = resp_div(dir_deriv, sign = T))



RDiv <- tibble(Div_loc_dir %>% dplyr::group_by(E1_ref, E2_ref) %>% 
                 summarise(mean = mean(div)))



m_RD <- gam(mean ~ te(E1_ref, E2_ref, k = c(5, 5)),
            data = RDiv,
            method = "REML")


RDiv$predicted_mean <- predict(m_RD, newdata = RDiv, type = "response")
# Generate grid of E1_ref and E2_ref values
e1_values <- seq(min(RDiv$E1_ref), max(RDiv$E1_ref), length.out = 100)
e2_values <- seq(min(RDiv$E2_ref), max(RDiv$E2_ref), length.out = 100)
grid <- expand.grid(E1_ref = e1_values, E2_ref = e2_values)

# Predict mean values on the grid
grid$pred_mean <- predict(m_RD, newdata = grid, type = "response")
library(ggsci)
# Plot predicted surface with colored background and black contour lines
capacity_divergence1 <- ggplot(grid, aes(x = E1_ref, y = E2_ref, z = pred_mean, fill = pred_mean)) +
  geom_tile(aes(fill = pred_mean)) +
  geom_contour(color = "black", linetype = "solid",  n = 10) +
  scale_fill_viridis(option = "inferno", , limits = c(0, 0.85), alpha = 0.85) +  # Use cm.colors palette with 100 colors
  labs(tag = "(a)",
       x = "E1", y = "E2", title = "Community 1", fill = "") +
  theme_classic() +
  theme(legend.position = "none") +
  theme(
    text = element_text(size = 30),  # Set the size of all text
    title = element_text(size = 30),  # Set the size of titles
    axis.title = element_text(size = 25),  # Set the size of axis titles
    axis.text = element_text(size = 25),  # Set the size of axis labels
    axis.ticks = element_line(size = 0.5),  # Set the size of axis ticks
    legend.title = element_text(size = 30),  # Set the size of legend title
    legend.text = element_text(size = 20)) 

### Community 2
#### Create community
E1_series <- seq(0, 50, 1)
E2_series <- seq(0, 50, 1)

## define new_data for following GAM fitting

new_data <- tibble(E1 =  seq(0, 50, 0.2),
                   E2 = seq(0, 50, 0.2))
## Set species parameters
species_pars <- list(a1_mean = 1e-3,
                     b1_mean = 0.02,
                     z1_mean = 20,
                     w1_mean = 10,
                     a2_mean = 1e-3,
                     b2_mean = 0.002,
                     z2_mean = 10,
                     w2_mean = 5,
                     zint_mean = 0,
                     sd_rate_mean = 0,
                     z1_range = NA,
                     z2_range = NA)


species_pars$z1_range <- 40
species_pars$z2_range <- 5

# create reproducible community
s <- 10
set.seed(2456)
comm <- get_community(s, species_pars)

# Visualise spp performances 
performance_comm3 <- multi_performance_plot(comm)
performance_comm3


refs <- crossing(E1 = seq(0, 50, 5),
                 E2 = seq(0, 50, 5))

nested_gams <- comm %>% 
  nest(cols =-species) %>% 
  mutate(
    gams = map(cols, ~ gam(rate ~  te(E1, E2),
                           data = .x,
                           method = "REML")),
    predicted = map(gams, ~ predict(.x, newdata = new_data))
  )

#list of gams
m_list <- (nested_gams$gams)
#list of spp names
my_spp_names <- (nested_gams$species)
# get partial derivatives
pd_list <- modify_depth(m_list, 1, ~ get_partials(., refs))

# from list to tibble
pd_spp <- tibble(
  E1_ref = map(pd_list, "E1"),
  E2_ref = map(pd_list, "E2"),
  pd_E1 = map(pd_list, "pd_E1"),
  pd_E2 = map(pd_list, "pd_E2")) %>%
  dplyr::mutate(sp = my_spp_names) %>%
  relocate(sp, E1_ref, E2_ref, pd_E1, pd_E2) %>%
  unnest(E1_ref, E2_ref, pd_E1, pd_E2)


radius <- 1
num_arrows <- 1000
pd_spp <- crossing(angle = rep(seq(0, 2*pi, length = num_arrows)),
                   E1_ref = pd_spp$E1_ref,
                   E2_ref =pd_spp$E2_ref) %>%
  full_join(pd_spp) %>%
  mutate(E1 = cos(angle) * radius,
         E2 = sin(angle) * radius,
         dir_deriv = E1 * pd_E1 + E2 * pd_E2,
         unit_vec_mag = sqrt(E1^2 + E2^2))


# filter useless stuff
pd_spp <- pd_spp %>% 
  dplyr::select(angle, sp, E1_ref, E2_ref, dir_deriv)

#
# Calculating diversity for each direction and location this should be always the first step
# Dissimilarity
Div_loc_dir <- pd_spp %>% 
  dplyr::group_by(E1_ref, E2_ref, angle) %>% 
  summarise(div = resp_div(dir_deriv, sign = T))



RDiv <- tibble(Div_loc_dir %>% dplyr::group_by(E1_ref, E2_ref) %>% 
                 summarise(mean = mean(div)))

m_RD <- gam(mean ~ te(E1_ref, E2_ref, k = c(5, 5)),
            data = RDiv,
            method = "REML")

RDiv$predicted_mean <- predict(m_RD, newdata = RDiv, type = "response")
# Generate grid of E1_ref and E2_ref values
e1_values <- seq(min(RDiv$E1_ref), max(RDiv$E1_ref), length.out = 100)
e2_values <- seq(min(RDiv$E2_ref), max(RDiv$E2_ref), length.out = 100)
grid <- expand.grid(E1_ref = e1_values, E2_ref = e2_values)

# Predict mean values on the grid
grid$pred_mean <- predict(m_RD, newdata = grid, type = "response")

# Plot predicted surface with colored background and black contour lines
capacity_divergence2 <- ggplot(grid, aes(x = E1_ref, y = E2_ref, z = pred_mean, fill = pred_mean)) +
  geom_tile(aes(fill = pred_mean)) +
  geom_contour(color = "black", linetype = "solid",  n = 10) +
  scale_fill_viridis(option = "inferno", , limits = c(0, 0.55), alpha = 0.85) +
  labs(tag = "(b)", title = "Community 2",
       x = "E1", y = "E2", fill = "Response\nCapacity") +
  theme_classic() +
  theme(
    text = element_text(size = 30),  # Set the size of all text
    title = element_text(size = 30),  # Set the size of titles
    axis.title = element_text(size = 25),  # Set the size of axis titles
    axis.text = element_text(size = 25),  # Set the size of axis labels
    axis.ticks = element_line(size = 0.5),  # Set the size of axis ticks
    legend.title = element_text(size = 30),  # Set the size of legend title
    legend.text = element_text(size = 20)) 




combined_plots <- capacity_divergence1 + capacity_divergence2

ggsave(combined_plots, file = here("Second version/figures_v2/Fig_3.png"), width = 20, height = 10)
ggsave(combined_plots, file = here("Second version/figures_v2/Fig_3.pdf"), width = 20, height = 10)



########## Figure 4 ##################
### first panel 

E1_series <- 273.15 + seq(0, 50, 1)
E2_series <- seq(0, 50, 1)


## Set parameter values - without interaction
pars <- list(a1 = 1e-9,
             b1 = 0.063,
             z1 = 285,
             w1 = 60,
             a2 = 1e-3,
             b2 = 0.02,
             z2 = 20,
             w2 = 10,
             z_int21 = 0,
             wint = 0,
             sd_rate = 0)
## Make a surface / experiment 
expt <- Make_expt(E1_series, E2_series, pars)
expt[[1]]$E1 <- expt[[1]]$E1 - 273.15

## Visualise the performance surface of a species
pc_figs <- Plot_performance_curves(expt[[1]])
fig1 <- pc_figs[[1]]+ theme_classic() + pc_figs[[2]] + theme_classic()
fig1 


pp1 <- fig1[[1]] + labs(x = "E1", y = "Growth Rate",  tag = "(a)") +
  theme_classic(base_size = 15) + 
  guides(color=guide_legend(title="E2")) +
  fig1[[2]] +
  theme_classic(base_size = 15) +
  labs(x = "E2", y = "Growth Rate", tag = "")+
  guides(color=guide_legend(title="E1"))

pp1 




E1_series <- seq(0, 50, 1)
E2_series <- seq(0, 50, 1)


## Set parameter values - without interaction
pars <- list(a1 = 1e-3,
             b1 = 0.02,
             z1 = 20,
             w1 = 10,
             a2 = 1e-3,
             b2 = 0.02,
             z2 = 20,
             w2 = 10,
             z_int21 = 0,
             wint = 0,
             sd_rate = 0)
## Make a surface / experiment 
expt <- Make_expt(E1_series, E2_series, pars)


## Visualise the performance surface of a species
pc_figs <- Plot_performance_curves(expt[[1]])
fig1 <- pc_figs[[1]]+ theme_classic() + pc_figs[[2]] + theme_classic()
fig1 



## Set parameter values - with interaction
pars <- list(a1 = 1e-3,
             b1 = 0.02,
             z1 = 20,
             w1 = 10,
             a2 = 1e-3,
             b2 = 0.02,
             z2 = 20,
             w2 = 10,
             z_int21 = 0.2,
             wint = 0,
             sd_rate = 0)
## Make a surface / experiment 
expt <- Make_expt(E1_series, E2_series, pars)


## Visualise the performance surface of a species
pc_figs <- Plot_performance_curves(expt[[1]])
fig2 <- pc_figs[[1]]+ theme_classic() + pc_figs[[2]] + theme_classic()

fig1 + theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5), 
        plot.title.position = "plot")
pp2 <- fig1[[1]] + labs(x = "E1", y = "Growth Rate") +
  theme_classic(base_size = 15) +  
  guides(color=guide_legend(title="E2")) +
  fig1[[2]]  +
  theme_classic(base_size = 15) +
  labs(x = "E2", y = "Growth Rate", tag = "")+
  guides(color=guide_legend(title="E1"))




pp1 <- pp1 + ggtitle("E1 dominant environmental variable") +
  theme(plot.title = element_text(hjust = 3))

pp2 <- pp2 + ggtitle("Equal effect of E1 and E2") +
  theme(plot.title = element_text(hjust = 100))

simulation_curves <- pp1/pp2 
simulation_curves
# create data frame
library(truncnorm)

set.seed(46354)
df_low_z1 <- data.frame(E1 = c(rtruncnorm(n=1000, a=0, b=50, mean=25, sd=3), rtruncnorm(n=1000, a=0, b=50, mean=24, sd=3), rtruncnorm(n=1000, a=0, b=50, mean=26, sd=3),
                               rtruncnorm(n=1000, a=0, b=50, mean=28, sd=3)),
                        species = rep(c("Sp 1", "Sp 2", "Sp 3", "Sp 4"), each = 1000))

# plot normal distributions
low_z1_div <- ggplot(df_low_z1, aes(x = E1, fill = species)) +
  annotate('rect', xmin=0, xmax=15, ymin=0.0, ymax=0.2, alpha=.2, fill='red')+
  annotate('text', x=8, y=0.205, label='Low mean E1')+
  annotate('rect', xmin=15, xmax=35, ymin=0.0, ymax=0.2, alpha=.2, fill='green')+
  annotate('text', x=23.5, y=0.205, label='Intermediate mean E1')+
  annotate('rect', xmin=35, xmax=50, ymin=0.0, ymax=0.2, alpha=.2, fill='blue')+
  annotate('text', x=42, y=0.205, label='High mean E1')+
  geom_density(alpha = 0.8) +
  scale_fill_viridis_d(option = "magma") +
  theme_classic() +
  labs(tag = "(c)", y = "Growth rate") + ggtitle("Low interspecific diversity in E1 optima") +labs(x = "E1") +
  theme(
    plot.title = element_text(hjust = 0.5)) +
  theme_classic(base_size = 15) +
  theme(
    plot.title = element_text(hjust = 0.5)) +
  xlim(0, 50) +
  theme(
    plot.title = element_text(hjust = 0.5)) + 
  theme(legend.position = "none")



# create data frame
set.seed(46354)
df_medium_z1 <- data.frame(E1 = c(rtruncnorm(n=1000, a=0, b=50, mean=15, sd=3), rtruncnorm(n=1000, a=0, b=50, mean=24, sd=3), rtruncnorm(n=1000, a=0, b=50, mean=30, sd=3),
                                  rtruncnorm(n=1000, a=0, b=50, mean=35, sd=3)),
                           species = rep(c("Sp 1", "Sp 2", "Sp 3", "Sp 4"), each = 1000))
# plot normal distributions
medium_z1_div <- ggplot(df_medium_z1, aes(x = E1, fill = species)) +
  annotate('rect', xmin=0, xmax=15, ymin=0.0, ymax=0.2, alpha=.2, fill='red')+
  annotate('text', x=8, y=0.205, label='Low mean E1')+
  annotate('rect', xmin=15, xmax=35, ymin=0.0, ymax=0.2, alpha=.2, fill='green')+
  annotate('text', x=23.5, y=0.205, label='Intermediate mean E1')+
  annotate('rect', xmin=35, xmax=50, ymin=0.0, ymax=0.2, alpha=.2, fill='blue')+
  annotate('text', x=42, y=0.205, label='High mean E1')+
  geom_density(alpha = 0.8) +
  scale_fill_viridis_d(option = "magma") +
  theme_classic() +
  labs(y = "Growth rate") + ggtitle("Medium interspecific diversity in E1 optima") +labs(x = "E1") +
  theme(
    plot.title = element_text(hjust = 0.5))+
  theme_classic(base_size = 15) +
  theme(
    plot.title = element_text(hjust = 0.5)) +
  xlim(0, 50) +
  theme(
    plot.title = element_text(hjust = 0.5)) + 
  theme(legend.position = "none")


medium_z1_div

# create data frame
set.seed(46354)
df_high_z1 <- data.frame(E1 = c(rtruncnorm(n=1000, a=0, b=50, mean=5, sd=4), rtruncnorm(n=1000, a=0, b=50, mean=18, sd=3), rtruncnorm(n=1000, a=0, b=50, mean=26, sd=3),
                                rtruncnorm(n=1000, a=0, b=50, mean=45, sd=3)),
                         species = rep(c("Sp 1", "Sp 2", "Sp 3", "Sp 4"), each = 1000))
                   
# plot normal distributions
high_z1_div <- ggplot(df_high_z1, aes(x = E1, fill = species)) +
  annotate('rect', xmin=0, xmax=15, ymin=0.0, ymax=0.2, alpha=.2, fill='red')+
  annotate('text', x=8, y=0.205, label='Low mean E1')+
  annotate('rect', xmin=15, xmax=35, ymin=0.0, ymax=0.2, alpha=.2, fill='green')+
  annotate('text', x=23.5, y=0.205, label='Intermediate mean E1')+
  annotate('rect', xmin=35, xmax=50, ymin=0.0, ymax=0.2, alpha=.2, fill='blue')+
  annotate('text', x=42, y=0.205, label='High mean E1')+
  geom_density(alpha = 0.8) +
  theme_classic(base_size = 15) +
  scale_fill_viridis_d(option = "magma") +
  theme_classic(base_size = 15) +
  labs(y = "Growth rate") + ggtitle("High interspecific diversity in E1 optima") +labs(x = "E1") +
  theme(
    plot.title = element_text(hjust = 0.5))+
  xlim(0, 50) +
  theme(
    plot.title = element_text(hjust = 0.5))

high_z1_div




comb <- low_z1_div / medium_z1_div / high_z1_div +
  theme(legend.position = "bottom")  # Set legend position for each plot

fig_optimum <- comb + 
  plot_layout( heights = c(0.5, 0.5, 0.6))
# Position of the optimum


scenario1 <- optima_scenario1 %>%  ggplot(aes(x = E1, y = E2, col = as.factor(community))) +
  geom_point(size = 8) +
  scale_colour_viridis_d()+
  labs(tag = "(b)", x = "Optimum E1", y = "Optimum E2",  col = ("community"), title = "Communities vary in only the amount of\n interspecific diversity in E1 optima \n - No correlation") +
  theme_classic(base_size = 15) +
  ylim(0, 50) +
  theme(axis.text = element_blank(),  # Remove axis labels
        axis.ticks = element_blank())


scenario2 <- optima_scenatio2 %>%  ggplot(aes(x = E1, y = E2, col = as.factor(community))) +
  geom_point(size = 8) +
  scale_colour_viridis_d()+
  labs(x = "Optimum E1", y = "Optimum E2", col = ("community"), title = "Communities vary in both the amount of\n interspecific diversity in E1 and E2 optima \n - Positive correlation") +
  theme_classic(base_size = 15) +
  ylim(0, 50) +
  theme(axis.text = element_blank(),  # Remove axis labels
        axis.ticks = element_blank())

scenario3 <- optima_scenario3 %>%  ggplot(aes(x = E1, y = E2, col = as.factor(community))) +
  geom_point(size = 8) +
  scale_colour_viridis_d()+
  labs(x = "Optimum E1", y = "Optimum E2", col = ("community"), title = "Communities vary in both the amount of\n interspecific diversity in E1 and E2 optima \n - Negative correlation") +
  theme_classic(base_size = 15) +
  ylim(0, 50) +
  theme(axis.text = element_blank(),  # Remove axis labels
        axis.ticks = element_blank())




optimum_scenarios <- scenario1 / scenario2 / scenario3 & theme(legend.position = "bottom")
optimum_scenarios <- optimum_scenarios + patchwork::plot_layout(guides = "collect")
optimum_scenarios



pp1 <- pp1 + ggtitle("E1 dominant environmental variable") +
  theme(plot.title = element_text(hjust = 2))

pp2 <- pp2 + ggtitle("Equal effect of E1 and E2") +
  theme(plot.title = element_text(hjust = 4.5))

simulation_curves <- pp1/pp2 
simulation_curves

left_panel <- (simulation_curves + fig_optimum) + plot_layout(heights = c(0.1, 0.1, 0.6))

# fig_4 <- left_panel + optimum_scenarios + plot_layout(ncol = 2)


fig_4 <- ggpubr::ggarrange(left_panel, optimum_scenarios, widths = c(3,2), 
                  ncol = 2)  
fig_4
# fig_Optimum_complete <- ggpubr::ggarrange(fig_optimum, optimum_scenarios, final_plot_corr, heights = c(2,2,2), 
#                                           ncol = 3)  
# fig_Optimum_complete 


ggsave(fig_4, file = here("Second version/figures_v2/Fig4.jpg"), width = 15, height = 18)
ggsave(fig_4, file = here("Second version/figures_v2/Fig4.png"), width = 15, height = 18)
ggsave(fig_4, file = here("Second version/figures_v2/Fig4.pdf"), width = 15, height = 18)




############ Figure 5 ################

# Figure 5
theme_set(theme_classic(base_size = 15))
div1 <-  ggplot(dd_divergence_scenario1, aes(x = z1, y = sign, color = E1_mean)) +
  # scale_y_continuous(limits = c(0.99, 1.005), expand = c(0.02, 0.02)) +
  scale_color_uchicago() +
  labs(x = "E1 optimum diversity", y = "Divergence (sign sensitive)") +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 16),
    panel.grid = element_blank()
  ) +
  geom_jitter(size = 2, alpha = 0.40, width = 0.2) +
  stat_summary(fun = mean, geom = "point", size = 5) + labs(tag = "(a)")+ 
  scale_x_discrete(limits = level_order)  + ggtitle("Only variation in E1\n optimum diversity")


div2 <-  ggplot(dd_divergence_scenario2, aes(x = z1, y = sign, color = E1_mean)) +
  # scale_y_continuous(limits = c(0.99, 1.005), expand = c(0.02, 0.02)) +
  scale_color_uchicago() +
  labs(x = "E1 optimum diversity", y = "Divergence (sign sensitive)") +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 16),
    panel.grid = element_blank()
  ) +
  geom_jitter(size = 2, alpha = 0.40, width = 0.2) +
  stat_summary(fun = mean, geom = "point", size = 5) + labs(tag = "(b)")+ 
  scale_x_discrete(limits = level_order) + ggtitle("Variation in E1 and E2\noptimum diversity - positive correlation")

div3 <-  ggplot(dd_divergence_scenario3, aes(x = z1, y = sign, color = E1_mean)) +
  # scale_y_continuous(limits = c(0.99, 1.005), expand = c(0.02, 0.02)) +
  scale_color_uchicago() +
  labs(x = "E1 optimum diversity", y = "Divergence (sign sensitive)") +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 16),
    panel.grid = element_blank()
  ) +
  geom_jitter(size = 2, alpha = 0.40, width = 0.2) +
  stat_summary(fun = mean, geom = "point", size = 5) + labs(tag = "(c)")+ 
  scale_x_discrete(limits = level_order) + ggtitle("Variation in E1 and E2\noptimum diversity - negative correlation")


div4 <-  ggplot(dd_divergence_scenario1_SAME, aes(x = z1, y = sign, color = E1_mean)) +
  # scale_y_continuous(limits = c(0.99, 1.005), expand = c(0.02, 0.02)) +
  scale_color_uchicago() +
  labs(x = "E1 optimum diversity", y = "Divergence (sign sensitive)") +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 16),
    panel.grid = element_blank()
  ) +
  geom_jitter(size = 2, alpha = 0.40, width = 0.2) +
  stat_summary(fun = mean, geom = "point", size = 5) + labs(tag = "(d)")+ 
  scale_x_discrete(limits = level_order)   + ggtitle("Only variation in E1\n optimum diversity")


div5 <-  ggplot(dd_divergence_scenario2_SAME, aes(x = z1, y = sign, color = E1_mean)) +
  # scale_y_continuous(limits = c(0.99, 1.005), expand = c(0.02, 0.02)) +
  scale_color_uchicago() +
  labs(x = "E1 optimum diversity", y = "Divergence (sign sensitive)") +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 16),
    panel.grid = element_blank()
  ) +
  geom_jitter(size = 2, alpha = 0.40, width = 0.2) +
  stat_summary(fun = mean, geom = "point", size = 5) + labs(tag = "(e)")+ 
  scale_x_discrete(limits = level_order) + ggtitle("Variation in E1 and E2\noptimum diversity - positive correlation")


div6 <-  ggplot(dd_divergence_scenario3_SAME, aes(x = z1, y = sign, color = E1_mean)) +
  # scale_y_continuous(limits = c(0.99, 1.005), expand = c(0.02, 0.02)) +
  scale_color_uchicago() +
  labs(x = "E1 optimum diversity", y = "Divergence (sign sensitive)") +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 16),
    panel.grid = element_blank()
  ) +
  geom_jitter(size = 2, alpha = 0.40, width = 0.2) +
  stat_summary(fun = mean, geom = "point", size = 5) + labs(tag = "(f)")+ 
  scale_x_discrete(limits = level_order) + ggtitle("Variation in E1 and E2\noptimum diversity - negative correlation")





top_plot      <- wrap_elements((div1 + div2 + div3 & theme(legend.position = "bottom") & 
                                  labs(col = "E1 mean value")) + plot_layout(guides = "collect")+
                                 plot_annotation(title = "E1 dominant environmental variable",
                                                 theme = theme(plot.title = element_text(hjust = 0.5))))


bottom_plot   <- wrap_elements((div4 + div5 + div6 & theme(legend.position = "bottom") & 
                                  labs(col = "E1 mean value")) + plot_layout(guides = "collect")+ 
                                 plot_annotation(title = "Equal effect of E1 and E2",
                                                 theme = theme(plot.title = element_text(hjust = 0.5))))

combined_final <- (top_plot / bottom_plot) 

combined_final <-combined_final + plot_layout(guides = "collect") # doesn't actually collect the legend... if you happen to know how to meke it work, please get in touch
combined_final <- combined_final
  theme(plot.caption = element_text(hjust = 0, vjust = 0.3, color = "blue", size = 15))


ggsave(combined_final, file = here("Second version/figures_v2/Fig5.jpg"), width = 18, height = 14)
ggsave(combined_final, file = here("Second version/figures_v2/Fig5.png"), width = 18, height = 14)
ggsave(combined_final, file = here("Second version/figures_v2/Fig5.pdf"), width = 18, height = 14)


