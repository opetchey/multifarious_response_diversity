# Manuscript figure

# rm(list = ls())

# packages
pkgs <- c("here", "dplyr", "ggplot2", "mgcv", "gratia", "scales", "tidyverse", "purrr")
vapply(pkgs, library, logical(1L), logical.return = TRUE, character.only = TRUE)

source.files <- list.files(here("r"), full.names = TRUE)
sapply(source.files, source, .GlobalEnv)

load(here("Data/my_work_space.RData"))
# the above loads all these old gratia functions - you should use the dev
# version of gratia and install the binary from my r-universe
# Now I delete all the old functions that are stored in the sourced image
rm(list = ls(getNamespace("gratia"), all.names = TRUE))
# this will throw warnings if you have a newer versio of gratia installed than
# is represented by the dump of R functions into the `./r` folder. This is
# harmless.

############ Figure 1 ################
#### Fig 1 normal distribution with two tangent lines ### 

# Set the mean and standard deviation
mu <- 50
sigma <- 25

# Create a data frame with x and y values for the normal distribution
df <- data.frame(x = seq(0, 100, length.out = 100))
df$y <- dnorm(df$x, mean = mu, sd = sigma)

# Plot the normal distribution using ggplot
p <- ggplot(df, aes(x, y)) +
  geom_line() 

# Add the first tangent line
x1 <- 25
y1 <- dnorm(x1, mean = mu, sd = sigma)
slope1 <- (dnorm(x1 + 0.01, mean = mu, sd = sigma) - dnorm(x1 - 0.01, mean = mu, sd = sigma)) / 0.02
intercept1 <- y1 - slope1 * x1
p <- p + geom_abline(intercept = intercept1, slope = slope1, col = "red")

# Add the second tangent line
x2 <- 75
y2 <- dnorm(x2, mean = mu, sd = sigma)
slope2 <- (dnorm(x2 + 0.01, mean = mu, sd = sigma) - dnorm(x2 - 0.01, mean = mu, sd = sigma)) / 0.02
intercept2 <- y2 - slope2 * x2
p <- p + geom_abline(intercept = intercept2, slope = slope2, col = "blue") +
  theme_classic(base_size = 14)  + theme(axis.title=element_text(size=15),
                          axis.text=element_text(size=12))+
  labs(x = "Environmental driver", y = "Growth rate", tag = "(a)")

# Display panel (a)
print(p)

# Creating panel (b)
m <- gam(y ~ s(x), data = df, method = "REML")

new_data <- data.frame(x = seq(0, 100, length.out = 1000))

# get first derivatives
d1_m <- gratia::derivatives(m, data = new_data)

d_plot <- ggplot(data = d1_m, mapping = aes(x = data, y = derivative)) +
  theme_classic(base_size = 14) + 
  labs(x = "Environmental driver",y = "Derivative",tag = "(b)") + 
  geom_hline(yintercept = 0,
             lty=2) +
  geom_line()

Fig1 <- p/d_plot
Fig1

# Save the plot
# ggsave(Fig1, file="/Users/francesco/Documents/GitHub/multifarious_response_diversity/Fig1.jpg", width=8, height=8)
# ggsave(Fig1, file="/Users/francesco/Documents/GitHub/multifarious_response_diversity/Fig1.pdf", width=8, height=8)
# Note you don't want JPEGs for images like this - that format is for pictures.
# You want PNG.
# Note I also modified the folder so these are dumped into the ms_figures folder
ggsave(Fig1, file = here("ms_figures/Fig1.jpg"), width = 8, height = 8)
ggsave(Fig1, file = here("ms_figures/Fig1.png"), width = 8, height = 8)
ggsave(Fig1, file = here("ms_figures/Fig1.pdf"), width = 8, height = 8)


############ Figure 2 ################
# Fig 2 - annotations added in power point



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

# Species 1 
df <- data_sim("eg2", n = 2000, dist = "normal", scale = 0.5, seed = 55)
# fit the GAM (note: for execution time reasons, k is set articifially low)
m <- gam(y ~  te(x, z, k = c(5, 5)), data = df, method = "REML")
# draw te(x,z)
p1 <- draw(m, rug = FALSE,
           n_contour = 30,
           contour_col = "lightgrey")


refs <- data.frame(x = c(0.10, 0.55),
                 z = c(0.45, 0.55))

refs <- get_partials(m, refs) %>%
  dplyr::rename(x_ref = x, z_ref = z)



radius <- 0.1
num_arrows <- 4
dd1 <- tibble(angle = rep(seq(0, 2*pi, length = num_arrows), nrow(refs))) %>%
  mutate(x_ref = rep(refs$x_ref, each = num_arrows),
         z_ref = rep(refs$z_ref, each = num_arrows)) %>%
  full_join(refs) %>%
  mutate(x = cos(angle) * radius,
         z = sin(angle) * radius,
         dir_deriv = x * pd_x + z * pd_z,
         unit_vec_mag = sqrt(x^2 + z^2))

col <- ifelse(dd1$x_ref==0.46, "A", 
              ifelse(dd1$z_ref==0.262, "A", 
                     "B"  )) # all other values map to NA

dd1$letter <- col

# need to replace p1 with smooth_estimates data and steal code from the `draw()`
# method to get the nice colours

# gs_p1 <- smooth_estimates(m) |>
#   ggplot(aes(x = x, y = z)) +
#   geom_raster(aes(fill = .estimate)) +
#   geom_contour(mapping = aes(z = .estimate),
#     colour = "lightgrey", bins = 30, na.rm = TRUE) +
#   scale_fill_distiller(palette = "RdBu", type = "div") +
#   expand_limits(fill = c(-1, 1) * max(abs(gs_p1_df[[".estimate"]]),
#     na.rm = TRUE))

Fig2 <- p1 +
  geom_segment(data = dd1,
    aes(x = x_ref, y = z_ref,
      xend = x_ref + x, yend = z_ref + z,
      col = dir_deriv), alpha = 0.9,
    linewidth = 3) +
  scale_colour_gradient2(low = muted("blue"),
    mid = "white",
    high = muted("red"),
    midpoint = 0) +
  geom_point(data = dd1,
    aes(x = x_ref, y = z_ref,
      size = 2)) +
  annotate("text", x = 0.10, y = 0.60, label = "A", size = 6) +
  annotate("text", x = 0.60, y = 0.70, label = "B", size = 6) +
  labs(x = "Temperature (k)", y = "Salinity (ppt)",
    color = "Directional\nDerivative") +
  theme_classic(base_size = 20) +
  # here is the issue - you need guide_colourbar() not guide_legend()
  guides(fill = guide_colourbar(title="Growth\nRate"), size = "none")

Fig2


# notice though that the directional deriative lines are not the same length!
# this is because the aspect ratio of the plot is not 1
# this versions fix that
ggsave(Fig2 + coord_equal(), file = here("ms_figures/Fig2-eq.jpg"), width = 10, height = 8)
ggsave(Fig2 + coord_equal(), file = here("ms_figures/Fig2-eq.png"), width = 10, height = 8)
ggsave(Fig2 + coord_equal(), file = here("ms_figures/Fig2-eq.pdf"), width = 10, height = 8)

############ Figure 3 ################
# Fig 3 is generated in the document "Appendix1_principles and demos" at line 668




############ Figure 4 ################
# Fig 4 
source.files <- list.files(here("r"), full.names = TRUE)
sapply(source.files, source, .GlobalEnv)

load(here("Data/my_work_space.RData"))
# the above loads all these old gratia functions - you should use the dev
# version of gratia and install the binary from my r-universe
# Now I delete all the old functions that are stored in the sourced image
rm(list = ls(getNamespace("gratia"), all.names = TRUE))
# this will throw warnings if you have a newer versio of gratia installed than
# is represented by the dump of R functions into the `./r` folder. This is
# harmless.


## Fig 4 - first part of the code also used for fig 5

library(pracma)

E1_series <- 273.15 + seq(0, 50, 1)
E2_series <- seq(0, 50, 1)
# Simulate spp performance curves with the modified Eppley function with and without interactive effect.

# Without interaction
set.seed(2465)
s <- 4 ## number of species
a1_mean <- 1e-9
b1_mean <- 0.063
z1_mean <- 285
w1_mean <- 60
a2_mean <- 1e-3
b2_mean <- 0.02
z2_mean <- 20
w2_mean <- 10
zint_mean <- 0 ## no interaction
sd_rate_mean <- 0.02
z1_range <- 40
z2_range <- 20
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
new_data <- tibble(E1 = 273.15 + seq(0, 50, 0.2),
                   E2 = seq(0, 50, 0.2))

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
refs2 <- crossing(E1 = 273.15 + seq(0, 50, 10),
                  E2 = seq(0, 50, 10))

#community 1 
m_list_no <- (nested_gams_no_inter$gams)

my_spp_names_no <- (nested_gams_no_inter$species)


(pd_list_no <- modify_depth(m_list_no, 1, ~ get_partials(., refs2)))




# create a dataframe with E1 and E2 changing over time
set.seed(64)

library(MASS)
refs_ts <- mvrnorm(n = 10,        #  specifies the sample size
                   mu = c(mean(E1_series), mean(E2_series)), #specifies the mean values of each column
                   Sigma = matrix(c( 105, 0.0,
                                     0.0, 105), #  specifies the correlation matrix
                                  nrow = 2))


refs_ts <- refs_ts %>%  as_tibble(refs_ts) %>% 
  dplyr::rename(E1 = "V1", E2 = "V2")

refs_ts <- refs_ts%>% 
  dplyr::mutate(time = seq.int(nrow(refs_ts)))

p_E1 <- refs_ts %>% 
  ggplot(aes(x = time, y = E1)) + geom_line( size = 1.5) +
  geom_point(size = 3)+
  theme_classic(base_size = 20) +
  labs(y = "Temperature", x = "Time", tag = "(a)") + 
  scale_x_continuous(breaks=seq(0, 10, 1))

p_E2 <- refs_ts %>% 
  ggplot(aes(x = time, y = E2)) + geom_line( size = 1.5) +
  geom_point(size = 3)+
  theme_classic(base_size = 20) +
  labs(y = "Salinity", x = "Time",tag = "(b)") + 
  scale_x_continuous(breaks=seq(0, 10, 1))

# Plot environmental change over time
p_E1 + p_E2
indipendent_fluctuations <- p_E1 + p_E2



# Generating spline for E1
xout <- seq(min(refs_ts$time), max(refs_ts$time), length.out = 100)
spline_E1 <- pracma::interp1(x = refs_ts$time, y = refs_ts$E1, xi = xout, method = "cubic")

# Creating a data frame for interpolated values
interpolated_df_E1 <- data.frame(time = xout, E1 = spline_E1)

# Adding the spline and points to the plot
p_E1_spline <- p_E1 +
  geom_line(data = interpolated_df_E1, aes(x = time, y = E1), color = "#1f77b4", size = 1) +
  geom_point(data = interpolated_df_E1, aes(x = time, y = E1), color = "#1f77b4", size = 3)+
  geom_point(size = 4)+ geom_line( size = 1.5)

# Viewing the plot
p_E1_spline


# Generating spline for E2
xout <- seq(min(refs_ts$time), max(refs_ts$time), length.out = 100)
spline_E2 <- pracma::interp1(x = refs_ts$time, y = refs_ts$E2, xi = xout, method = "cubic")

# Creating a data frame for interpolated values
interpolated_df_E2 <- data.frame(time = xout, E2 = spline_E2)

# Adding the spline and points to the plot
p_E2_spline <- p_E2 +
  geom_line(data = interpolated_df_E2, aes(x = time, y = E2), color = "#ff7f0e", size = 1) +
  geom_point(data = interpolated_df_E2, aes(x = time, y = E2), color = "#ff7f0e", size = 3)+
  geom_point(size = 4)+ geom_line( size = 1.5)

# Viewing the plot
p_E2_spline


# creating interpolated df of env condition
interpolated_df <- cbind(interpolated_df_E1, interpolated_df_E2[,2]) %>% rename(E2 = "interpolated_df_E2[, 2]")

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


dir_plot <- red_spp %>% dplyr::filter(sp == "s1") %>% 
  ggplot(., mapping = aes(x = time, y = (dir_deriv))) +
  theme_bw(base_size = 20)+
  geom_line()+
  geom_point(size = 3)+
  labs(x = "time",y = "Directional derivative", tag = "(d)") +
  geom_hline(yintercept = 0, linetype= "dashed") + theme_classic(base_size = 15) +
  scale_x_continuous(breaks=seq(0, 10, 1))





#### Using interpolated data
# get partial derivatives - measured values
(pd_list_interpolated <- modify_depth(m_list_no, 1, ~ get_partials(., interpolated_df)))
# from list to tibble
(pd_spp_interpolated <- tibble(
  E1_ref = map(pd_list_interpolated, "E1"),
  E2_ref = map(pd_list_interpolated, "E2"),
  pd_E1 = map(pd_list_interpolated, "pd_E1"),
  pd_E2 = map(pd_list_interpolated, "pd_E2")) %>%
    dplyr::mutate(sp = my_spp_names_no) %>%
    relocate(sp, E1_ref, E2_ref, pd_E1, pd_E2) %>%
    unnest(E1_ref, E2_ref, pd_E1, pd_E2))
# add time
(pd_spp_interpolated <-cbind(pd_spp_interpolated, interpolated_df$time) %>% 
    dplyr::rename(time = "interpolated_df$time"))

# calculation next value for directional derivatives, and get directional derivatives
(pd_spp_interpolated <- pd_spp_interpolated %>% transform( nxt_value_E1 = c(E1_ref[-1], NA)) %>%
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
red_interpolated <- pd_spp_interpolated %>%  dplyr::select(sp, time, E1_ref, E2_ref, dir_deriv)


dir_plot_interpolated <- red_interpolated %>%
  dplyr::filter(sp == "s1") %>%
  ggplot(aes(x = time, y = dir_deriv)) +
  theme_bw(base_size = 20) +
  geom_line(color = "#009E73", position = position_nudge(x = 0.5), size = 1) +
  geom_point(color = "#009E73", position = position_nudge(x = 0.5), size = 3) +
  labs(x = "time", y = "Directional derivative", tag = "(d)") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_classic(base_size = 20) +
  scale_x_continuous(breaks = seq(0.5, 10.5, 1), labels = seq(0, 10, 1), expand = c(0, 0.5))


red_spp_s1 <- red_spp %>% dplyr::filter(sp == "s1") 


# Adding data from pd_spp_interpolated
dir_plot <- dir_plot_interpolated +
  geom_line(data = red_spp_s1, aes(x = time + 0.5, y = dir_deriv), color = "black", size = 1) +
  geom_point(data = red_spp_s1, aes(x = time + 0.5, y = dir_deriv), color = "black", size = 3) +
  scale_x_continuous(breaks = seq(0, 10, 1), expand = c(0, 0.5))


# View combined plot
dir_plot

# Sp response surface

sp_plot <- sp1 + geom_point(data = refs_ts, mapping = aes(x = E1, y = E2), size = 3) +
  geom_path(data = refs_ts, size = 1.5) +
  geom_point(data = interpolated_df, mapping = aes(x = E1, y = E2), size = 1) +
  geom_path(data = interpolated_df,  size = 0.5, alpha = 0.5) + guides(color="none" , fill=guide_colourbar(title="Growth \n Rate")) +
  geom_label(data = refs_ts, mapping = aes(label = time), size = 10) + labs(title = "", tag = "(c)")+ theme_classic(base_size = 20) + ylab("Salinity (ppt)") + xlab("Temperature (k)")  + 
  geom_path(data = interpolated_df, size = 0.5, alpha = 0.5) 




Fig4 <- (p_E1_spline + p_E2_spline + sp_plot) / dir_plot
Fig4


ggsave(Fig4, file = here("ms_figures/Fig4.jpg"), width = 22, height = 15)
ggsave(Fig4, file = here("ms_figures/Fig4.png"), width = 22, height = 15)
ggsave(Fig4, file = here("ms_figures/Fig4.pdf"), width = 22, height = 15)


############ Figure 5 ################
#Fig 5
okabe_ito_palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")


# community 1 
# reduce the dataframe and keep only what we need
red_spp <- pd_spp_no %>%  dplyr::select(sp, time, E1_ref, E2_ref, dir_deriv)
# from long to wide
rdiv_sim_nocor <- red_spp %>%
  spread( sp, dir_deriv)

rdiv_sim_nocor[is.na(rdiv_sim_nocor)] <- 0


# actual calculation for only the same species used above
rdiv_sim_nocor$rdiv<-apply(rdiv_sim_nocor[,c(4:7)], 1, resp_div, sign_sens = F)
rdiv_sim_nocor$sign<-apply(rdiv_sim_nocor[,c(4:7)], 1, resp_div, sign_sens = T)
rdiv_sim_nocor$Med<-median(rdiv_sim_nocor$rdiv)
rdiv_sim_nocor



# community 1 
rmax<-max(rdiv_sim_nocor$rdiv); rmin<-1
dmax<-max(rdiv_sim_nocor$sign); dmin<-0
# community 1 


dd_plot <- ggplot(data = red_spp, mapping = aes(x = time, y = (dir_deriv), col = sp)) +
  theme_bw(base_size = 12)+
  geom_line(size = 1.5)+
  labs(x = "time",y = "Directional derivative", tag = "a")+
  geom_hline(yintercept = 0, linetype= "dashed")+
  scale_x_continuous(breaks=seq(0, 15, 1)) +
  scale_color_manual(values = okabe_ito_palette)

rdiv_plot <- ggplot(data = rdiv_sim_nocor, mapping = aes(x = time, y = rdiv)) +
  theme_bw(base_size = 12)+
  geom_line(size = 1.5)+
  labs(x = "time",y = "Dissimilarity (derivatives)", tag = "c") +
  theme_bw(base_size = 12)+
  geom_richtext(x = 2,
                mapping = aes(y = rmax - 0.005),
                size=4.5,
                label.color = NA,
                label = paste("Mean =",paste0(round(mean(rdiv_sim_nocor$rdiv),digits = 2))))+
  scale_x_continuous(breaks=seq(0, 15, 1))
  



sing_plot <- ggplot(data = rdiv_sim_nocor, mapping = aes(x = time, y = sign)) +
  theme_bw(base_size = 12)+
  geom_line(size = 1.5)+
  labs(x = "time",y = "Divergence", tag = "e") +
  theme_bw(base_size = 12)+
  geom_richtext(x = 2,
                mapping = aes(y = dmax - 0.1),
                size=4.5,
                label.color = NA,
                label = paste("Mean =",paste0(round(mean(rdiv_sim_nocor$sign),digits = 2))))+
  lims(y = c(dmin,dmax))+
  scale_x_continuous(breaks=seq(0, 15, 1))

sp1_nocor <- sp1_nocor + theme_classic(base_size = 15) + labs(title = "Species 1") + xlab("Temperature (k)") + ylab("Salinity (ppt)") + guides(fill=guide_colourbar(title="Growth \n Rate"))
sp2_nocor<- sp2_nocor + theme_classic(base_size = 15) + labs(title = "Species 2")+ xlab("Temperature (k)") + ylab("Salinity (ppt)") + guides(fill=guide_colourbar(title="Growth \n Rate"))
sp3_nocor <- sp3_nocor + theme_classic(base_size = 15) + labs(title = "Species 3")+ xlab("Temperature (k)") + ylab("Salinity (ppt)") + guides(fill=guide_colourbar(title="Growth \n Rate"))
sp4_nocor <- sp4_nocor + theme_classic(base_size = 15) + labs(title = "Species 4")+ xlab("Temperature (k)") + ylab("Salinity (ppt)") + guides(fill=guide_colourbar(title="Growth \n Rate"))

dd_plot <- dd_plot + labs(tag = "(e)")+ theme_classic(base_size = 15) + theme(legend.position = "top")
rdiv_plot <- rdiv_plot + labs(title = "", tag = "(f)")+ theme_classic(base_size = 15) + ylab("Dissimilarity")
sing_plot <- sing_plot + labs(tag = "(g)")+ theme_classic(base_size = 15)

try <- ggarrange(sp1_nocor, sp2_nocor, sp3_nocor,  sp4_nocor,  nrow = 2, ncol=2, common.legend = TRUE, legend="right")

fig5 <- ggarrange(try, dd_plot, rdiv_plot, sing_plot, heights = c(3.5, 1.5, 1.2, 1), widths = c(c(3.5, 1.5, 1.5, 1.5)),
                  ncol = 1, nrow = 4, align = "v") 
fig5


ggsave(fig5, file = here("ms_figures/Fig5.jpg"), width = 8, height = 12)
ggsave(fig5, file = here("ms_figures/Fig5.png"), width = 8, height = 12)
ggsave(fig5, file = here("ms_figures/Fig5.pdf"), width = 8, height = 12)

############ Figure 6 ################
# Fig 6 is generated in the document called "Appendix1_principle and demos.Rmd", line 707





############ Figure 7 ################
### Figure 7 - Creating the basis for the figure. Additional text and lines have been added in power point

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

# Species 1 
df <- data_sim("eg2", n = 2000, dist = "normal", scale = 0.5, seed = 55)
# fit the GAM (note: for execution time reasons, k is set articifially low)
m <- gam(y ~ ti(x) + ti(z) + te(x, z, k = c(5, 5)), data = df, method = "REML")
# draw te(x,z)
p1 <- draw(m, rug = FALSE,
           n_contour = 30,
           contour_col = "lightgrey")
refs <- crossing(x = seq(0.1, 0.9, 0.35),
                 z = seq(0.1, 0.9, 0.35))
refs <- get_partials(m, refs) %>%
  dplyr::rename(x_ref = x, z_ref = z)
radius <- 0.1
num_arrows <- 4
dd1 <- tibble(angle = rep(seq(0, 2*pi, length = num_arrows), nrow(refs))) %>%
  mutate(x_ref = rep(refs$x_ref, each = num_arrows),
         z_ref = rep(refs$z_ref, each = num_arrows)) %>%
  full_join(refs) %>%
  mutate(x = cos(angle) * radius,
         z = sin(angle) * radius,
         dir_deriv = x * pd_x + z * pd_z,
         unit_vec_mag = sqrt(x^2 + z^2))
sp1 <- p1 +
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
  xlab("Temperature (k)") + ylab("Salinity (ppt)") + 
  labs(title = "Sp 1", tag = "(a)", caption = "") + theme_classic(base_size = 20) +
  guides(size = "none") 

# Species 2
df <- data_sim("eg2", n = 2000, dist = "normal", scale = 0.5, seed = 55)
# fit the GAM (note: for execution time reasons, k is set articifially low)
m <- gam(y ~  te(x, z, k = c(5, 5)), data = df, method = "REML")
# draw te(x,z)
p1 <- draw(m, rug = FALSE,
           n_contour = 30,
           contour_col = "lightgrey")
refs <- crossing(x = seq(0.1, 0.9, 0.35),
                 z = seq(0.1, 0.9, 0.35))
refs <- get_partials(m, refs) %>%
  dplyr::rename(x_ref = x, z_ref = z)
radius <- 0.1
num_arrows <- 4
dd1 <- tibble(angle = rep(seq(0, 2*pi, length = num_arrows), nrow(refs))) %>%
  mutate(x_ref = rep(refs$x_ref, each = num_arrows),
         z_ref = rep(refs$z_ref, each = num_arrows)) %>%
  full_join(refs) %>%
  mutate(x = cos(angle) * radius,
         z = sin(angle) * radius,
         dir_deriv = x * pd_x + z * pd_z,
         unit_vec_mag = sqrt(x^2 + z^2))
sp2 <- p1 + 
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
  xlab("Temperature (k)") + ylab("Salinity (ppt)") + 
  labs(title = "Sp 2", tag = "(b)", caption = "") + theme_classic(base_size = 20) +
  guides(size = "none")

Fig_7 <- ggarrange(sp1, sp2, heights = c(2,2), 
          ncol = 2) 
Fig_7

ggsave(Fig_7, file = here("ms_figures/Fig7.jpg"), width = 14, height = 8)
ggsave(Fig_7, file = here("ms_figures/Fig7.png"), width = 14, height = 8)
ggsave(Fig_7, file = here("ms_figures/Fig7.pdf"), width = 14, height = 8)

### Preparing Figure 8 - Temperature dominant env variable on species growth rate
## define the series of values of the environmental variables. (We also provie the code to visualise the spp responses when we insert an interaction)
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


## Visualise the performance surface of a species
pc_figs <- Plot_performance_curves(expt[[1]])
fig1 <- pc_figs[[1]]+ theme_classic() + pc_figs[[2]] + theme_classic()
fig1 



## Set parameter values - with interaction
pars <- list(a1 = 1e-9,
             b1 = 0.063,
             z1 = 285,
             w1 = 60,
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

fig1 + theme_classic()
pp1 <- fig1[[1]] + labs(x = "Temperature", y = "Growth Rate",  tag = "(a)") +
  theme_classic(base_size = 15) +
  guides(color=guide_legend(title="Salinity")) +
  fig1[[2]] +
  theme_classic(base_size = 15) +
  labs(x = "Salinity", y = "Growth Rate", tag = "")+
  guides(color=guide_legend(title="Temperature"))

pp1 <- pp1 +  ggtitle("Temperature dominant environmental variable")+
  theme(plot.title = element_text(hjust = 1.5))


pp2 <- fig2[[1]] +  labs(x = "Temperature", y = "Growth Rate", tag = "(b)") +
  theme_classic(base_size = 15) +
  guides(color=guide_legend(title="Salinity")) +
  fig2[[2]] +
  theme_classic(base_size = 15) +
  labs(x = "Salinity", y = "Growth Rate", tag = "")+
  guides(color=guide_legend(title="Temperature"))

pp2 <- pp2 +  ggtitle("Equal effect of temperature and salinity")+
  theme(plot.title = element_text(hjust = 1))

simulation_curves_Tdom <- pp1/pp2
simulation_curves_Tdom





### Preparing Figure 8 - same effect of temperature and salinity on species growth rate
## define the series of values of the environmental variables
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
expt[[1]]$E1 <- expt[[1]]$E1 + 273.15

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
expt[[1]]$E1 <- expt[[1]]$E1 + 273.15

## Visualise the performance surface of a species
pc_figs <- Plot_performance_curves(expt[[1]])
fig2 <- pc_figs[[1]]+ theme_classic() + pc_figs[[2]] + theme_classic()

fig1 + theme_classic() 
pp2 <- fig1[[1]] + labs(x = "Temperature", y = "Growth Rate",  tag = "(b)") +
  theme_classic(base_size = 15) +
  guides(color=guide_legend(title="Salinity")) +
  fig1[[2]]  +
  theme_classic(base_size = 15) +
  labs(x = "Salinity", y = "Growth Rate", tag = "")+
  guides(color=guide_legend(title="Temperature"))

pp2 <- pp2 +  ggtitle("Equal effect of temperature and salinity")+
  theme(plot.title = element_text(hjust = 2))


simulation_curves <- pp1/pp2
simulation_curves


ggsave(simulation_curves, file = here("ms_figures/Fig8.jpg"), width = 10, height = 8)
ggsave(simulation_curves, file = here("ms_figures/Fig8.png"), width = 10, height = 8)
ggsave(simulation_curves, file = here("ms_figures/Fig8.pdf"), width = 10, height = 8)




############ Figure 9 ################
### Preparing Figure 9
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
  labs(y = "Growth rate") + ggtitle("Low interspecific diversity \n in temperature optima") +labs(x = "Temperature", tag = "(a)") +
  theme(
    plot.title = element_text(hjust = 0.5)) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5))



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
  labs(y = "Growth rate") + ggtitle("Medium interspecific diversity \n in temperature optima") +labs(x = "Temperature", tag = "(d)") +
  theme(
    plot.title = element_text(hjust = 0.5))+
  theme_classic(base_size = 12) +
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
  labs(y = "Growth rate") + ggtitle("High interspecific diversity \n in temperature optima") +labs(x = "Temperature", tag = "(g)") +
  theme(
    plot.title = element_text(hjust = 0.5))+
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5))
high_z1_div

comb <- low_z1_div/ medium_z1_div / high_z1_div & theme(legend.position = "bottom")
fig_optimum <- comb + patchwork::plot_layout(guides = "collect", heights = c(0.5, 0.5, 0.5))

# Position of the optimum
scenario1 <- optima_scenario1 %>%  ggplot(aes(x = E1, y = E2, col = as.factor(community))) +
  geom_point(size = 8) +
  scale_colour_viridis_d()+
  labs(x = "Optimum Temperature", y = "Optimum Salinity", tag = '(b)', col = ("community"), title = "Communities vary in only \n the amount of interspecific diversity \n in temperature optima
") +
  theme_classic(base_size = 12) + 
  ylim(0, 50) 



scenario2 <- optima_scenatio2 %>%  ggplot(aes(x = E1, y = E2, col = as.factor(community))) +
  geom_point(size = 8) +
  scale_colour_viridis_d()+
  labs(x = "Optimum Temperature", y = "Optimum Salinity", tag = '(e)',  col = ("community"), title = "Communities vary in both \n the amount of interspecific diversity \n in temperature and salinity optima \n - Positive correlation") +
  theme_classic(base_size = 12) + 
  ylim(0, 50) 

scenario3 <- optima_scenario3 %>%  ggplot(aes(x = E1, y = E2, col = as.factor(community))) +
  geom_point(size = 8) +
  scale_colour_viridis_d()+
  labs(x = "Optimum Temperature", y = "Optimum Salinity", tag = '(h)',  col = ("community"), title = "Communities vary in both \n the amount of interspecific diversity \n in temperature and salinity optima \n - Negative correlation") +
  theme_classic(12)+ 
  ylim(0, 50) 


# Generating simulated data for three scenarios
set.seed(123)

# Function to create a scatter plot with a regression line and confidence intervals
plot_scenario <- function(data) {
  ggplot(data, aes(x = Salinity, y = Temperature)) +
    geom_point() +
    xlim(0, 50) + ylim(293, 330) +  # Set x and y axis limits
    labs(x = "Salinity", y = "Temperature") +
    geom_smooth(method = "lm", se = TRUE) +  # Add correlation line with CI
    theme_minimal()  # Optional: Adjust plot theme as needed
}

# Scenario 1: No correlation
scenario_1 <- data.frame(
  Salinity = runif(30, min = 20, max = 30),
  Temperature = runif(30, min = 293, max = 330)
)


# Scenario 2: Positive correlation with adjusted temperature range
Salinity_2 <- runif(30, min = 0, max = 50)
scenario_2 <- data.frame(
  Salinity = Salinity_2,
  Temperature = Salinity_2 * 0.7 + runif(30, min = 290, max = 310)
)


# Scenario 3: Negative correlation
Salinity_3 <- runif(30, min = 0, max = 50)
scenario_3 <- data.frame(
  Salinity = Salinity_3,
  Temperature = -Salinity_3 * 0.7 + runif(30, min = 320, max = 330)
)

# Create plots for each scenario
plot_1 <- plot_scenario(scenario_1) + labs(x = "Diversity in optima position for Temperature", y = "Diversity in optima position \n for Salinity",
                                           title = "No Correlation", tag = '(c)') + theme(
                                             plot.title = element_text(hjust = 0.5)) + xlim(20, 30) + theme_classic()
plot_2 <- plot_scenario(scenario_2) + labs(x = "Diversity in optima position for Temperature", y = "Diversity in optima position \n for Salinity",
                                          title = "Positive Correlation", tag = '(f)') + theme(
                                            plot.title = element_text(hjust = 0.5))+ theme_classic()
plot_3 <- plot_scenario(scenario_3) +  labs(x = "Diversity in optima position for Temperature", y = "Diversity in optima position \n for Salinity",
                                            title = "Negative Correlation", tag = '(i)') + theme(
                                              plot.title = element_text(hjust = 0.5))+  theme_classic()

# Arrange plots in a column
final_plot_corr <- plot_1 / plot_2 / plot_3
final_plot_corr



optimum_scenarios <- scenario1 / scenario2 / scenario3 & theme(legend.position = "bottom")
optimum_scenarios <- optimum_scenarios + patchwork::plot_layout(guides = "collect")
optimum_scenarios

fig_Optimum_complete <- ggpubr::ggarrange(fig_optimum, optimum_scenarios, final_plot_corr, heights = c(2,2,2), 
                                  ncol = 3, align = "hv")  
fig_Optimum_complete


ggsave(fig_Optimum_complete, file = here("ms_figures/Fig9.jpg"), width = 15, height = 15)
ggsave(fig_Optimum_complete, file = here("ms_figures/Fig9.png"), width = 15, height = 15)
ggsave(fig_Optimum_complete, file = here("ms_figures/Fig9.pdf"), width = 15, height = 15)


############ Figure 10 ################
### Preparing Figure 10
# Starting with time series of temperature changing over time with different mean values (low, medium, high)
## define the series of values of the environmental variables
E1_series <- 210.15 + seq(0, 120, 1)
E1min = 180
E1max = 350



# Scenario 1 - low mean vakues of E1 and E2
theme_set(theme_bw(base_size = 12))
# create random distribution
set.seed(65354)
refs1 <- as.data.frame(rnorm(50, quantile(E1_series, .01 ), 10))
refs1$time <- seq.int(nrow(refs1))
refs1 <- refs1 %>% dplyr::rename(E1 = "rnorm(50, quantile(E1_series, 0.01), 10)")
# plot env change over time
p_E1 <- refs1 %>% 
  ggplot(aes(x = time, y = E1)) + geom_line() +
  theme_classic(base_size = 12) +
  geom_hline(yintercept = mean(refs1$E1),lty=2)+ 
  labs(tag = "(a)", y = "Temperature") +
  lims(y = c(E1min, E1max))

p_E1


### scenario 2 - intermediate mean values E1
theme_set(theme_bw(base_size = 12))
# create random distribution
set.seed(65354)
refs2 <- as.data.frame(rnorm(50, quantile(E1_series, 0.50), 10))
refs2$time <- seq.int(nrow(refs1))
refs2 <- refs2 %>% dplyr::rename(E1 = "rnorm(50, quantile(E1_series, 0.5), 10)")
# plot env change over time
p2_E1 <- refs2 %>% 
  ggplot(aes(x = time, y = E1)) + geom_line() +
  theme_classic(base_size = 12) +
  geom_hline(yintercept = mean(refs2$E1),lty=2)+ 
  labs(tag = "(b)", y = "Temperature") +
  lims(y = c(E1min, E1max))

p2_E1


### scenario 3 - high mean values E1
theme_set(theme_bw(base_size = 12))
# create random distribution
set.seed(65354)
refs3 <- as.data.frame(rnorm(50, quantile(E1_series, .99 ), 10))
refs3$time <- seq.int(nrow(refs1))
refs3 <- refs3 %>% dplyr::rename(E1 = "rnorm(50, quantile(E1_series, 0.99), 10)")
# plot env change over time
p3_E1 <- refs3 %>% 
  ggplot(aes(x = time, y = E1)) + geom_line() +
  theme_classic(base_size = 12) +
  geom_hline(yintercept = mean(refs3$E1),lty=2)+ 
  labs(tag = "(c)", y = "Temperature") +
  lims(y = c(E1min, E1max))

p3_E1


env <- (p_E1 + p2_E1 + p3_E1)



# Now conceptual figure on diversity of spp responses
# create data frame
library(truncnorm)

set.seed(46354)
df_low_z1 <- data.frame(E1 = c(rtruncnorm(n=1000, a=250, b=350, mean=270, sd=15), rtruncnorm(n=1000, a=250, b=350, mean=295, sd=15), rtruncnorm(n=1000, a=250, b=350, mean=272, sd=15),
                               rtruncnorm(n=1000, a=50, b=350, mean=255, sd=15)),
                        species = rep(c("Sp 1", "Sp 2", "Sp 3", "Sp 4"), each = 1000))

# plot normal distributions
low_z1_div <- ggplot(df_low_z1, aes(x = E1, fill = species)) +
  geom_density(alpha = 0.5) +
  scale_fill_viridis_d() +
  theme_classic() +
  labs(y = "Growth rate") + ggtitle("Low interspecific diversity in temperature optima") +labs(x = "Temperature", tag = "(d)") +
  theme(
    plot.title = element_text(hjust = 0.5)) +
  theme_classic(base_size = 15) +
  annotate('rect', xmin=250, xmax=305, ymin=0.0, ymax=0.035, alpha=.2, fill='red')+
  annotate('text', x=280, y=0.032, label='Intermediate mean Temperature')+
  annotate('rect', xmin=200, xmax=250, ymin=0.0, ymax=0.035, alpha=.2, fill='green')+
  annotate('text', x=220, y=0.032, label='Low mean Temperature')+
  annotate('rect', xmin=305, xmax=350, ymin=0.0, ymax=0.035, alpha=.2, fill='blue')+
  annotate('text', x=330, y=0.032, label='High mean Temperature')+
  theme(
    plot.title = element_text(hjust = 0.5))




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
  labs(y = "Growth rate") + ggtitle("High interspecific diversity in temperature optima") +labs(x = "Temperature", tag = "(e)") +
  theme(
    plot.title = element_text(hjust = 0.5))+
  theme_classic(base_size = 15) +
  annotate('rect', xmin=250, xmax=305, ymin=0.0, ymax=0.035, alpha=.2, fill='red')+
  annotate('text', x=280, y=0.032, label='Intermediate mean Temperature')+
  annotate('rect', xmin=200, xmax=250, ymin=0.0, ymax=0.035, alpha=.2, fill='green')+
  annotate('text', x=220, y=0.032, label='Low mean Temperature')+
  annotate('rect', xmin=305, xmax=350, ymin=0.0, ymax=0.035, alpha=.2, fill='blue')+
  annotate('text', x=330, y=0.032, label='High mean Temperature') +
  theme(
    plot.title = element_text(hjust = 0.5))


comb <- low_z1_div / high_z1_div & theme(legend.position = "bottom")
fig_optimum <- comb + plot_layout(guides = "collect", heights = c(0.5, 0.5))

fig_optimum_color <- comb + plot_layout(guides = "collect")

mean_value_env_change <- ggarrange(env, fig_optimum_color, ncol = 1, nrow = 2, heights = c(0.3, 0.7), align = "v")

mean_value_env_change

ggsave(mean_value_env_change, file = here("ms_figures/Fig10.jpg"), width = 10, height = 9)
ggsave(mean_value_env_change, file = here("ms_figures/Fig10.png"), width = 10, height = 9)
ggsave(mean_value_env_change, file = here("ms_figures/Fig10.pdf"), width = 10, height = 9)


############ Figure 11 ################
# Figure 11 - potential response diversity
# Potential response diversity - final overview
col2 <- ifelse(Absolute_RD$scenario==1, "None",
               ifelse(Absolute_RD$scenario==2, "Positive",
                      ifelse(Absolute_RD$scenario==3, "Negative",
                                    NA  ))) # all other values map to NA

Absolute_RD$correlation <- col2
Absolute_RD$correlation <- factor(Absolute_RD$correlation, levels = c("None", "Positive", "Negative"))
Absolute_plot_tot <-  ggplot(Absolute_RD, aes(x = correlation, y = value, color = community)) +
  scale_color_uchicago(labels=c('High', 'Medium', 'Low')) + 
  labs(x = "Correlation in diversity of responses", y = "Response Diversity") + facet_wrap("RDiv",  scales = "free_y") +
   facet_wrap("RDiv",  scales = "free_y") +
  geom_jitter( size = 6, alpha = 0.80, width = 0.20) +
  theme_bw(base_size = 12)  + ggtitle("Treatment I:\nTemperature dominant environmental variable") +
  guides(color=guide_legend(title="Treatment II:\nDiversity in responses to temperature"))


col2 <- ifelse(Absolute_RD_same$scenario==1, "None",
               ifelse(Absolute_RD_same$scenario==2, "Positive",
                      ifelse(Absolute_RD_same$scenario==3, "Negative",
                             NA  ))) # all other values map to NA

Absolute_RD_same$correlation <- col2
Absolute_RD_same$correlation <- factor(Absolute_RD_same$correlation, levels = c("None", "Positive", "Negative"))

Absolute_plot_tot_same <-  ggplot(Absolute_RD_same, aes(x = correlation, y = value, color = community)) +
  scale_color_uchicago(labels=c('High', 'Medium', 'Low')) + 
  labs(x = "Correlation in diversity of responses", y = "Response Diversity") + facet_wrap("RDiv",  scales = "free_y") +
  facet_wrap("RDiv",  scales = "free_y") +
  geom_jitter( size = 6, alpha = 0.80, width = 0.20) +
  theme_bw(base_size = 12) +
  guides(color=guide_legend(title="Treatment II:\nDiversity in responses to temperature"))

Absolute_RD_ms <-( Absolute_plot_tot + labs(tag = "(a)"))  /( Absolute_plot_tot_same + labs(tags = "(b)")) & theme(legend.position = "right")
Absolute_RD_ms <- Absolute_RD_ms + plot_layout(guides = "collect") + ggtitle("Treatment I:\nEqual effect of temperature and salinity")

Absolute_RD_ms

ggsave(Absolute_RD_ms, file = here("ms_figures/Fig11.jpg"), width = 10, height = 10)
ggsave(Absolute_RD_ms, file = here("ms_figures/Fig11.png"), width = 10, height = 10)
ggsave(Absolute_RD_ms, file = here("ms_figures/Fig11.pdf"), width = 10, height = 10)


############ Figure 12 ################
# Figure 12
theme_set(theme_classic(base_size = 15))
diss1 <-  ggplot(dd_dissimilarity_scenario1, aes(x = z1, y = rdiv, color = E1_mean)) +
  # scale_y_continuous(limits = c(0.99, 1.005), expand = c(0.02, 0.02)) +
  scale_color_uchicago() +
  labs(x = "Temperature optimum diversity", y = "Dissimilarity (derivatives)") +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 16),
    panel.grid = element_blank()
  ) +
  geom_jitter(size = 2, alpha = 0.25, width = 0.2) +
  stat_summary(fun = mean, geom = "point", size = 5) + labs(tag = "(a)")+ 
  scale_x_discrete(limits = level_order)  + ggtitle("Treatment II:\nOnly variation in Temperature\n optimum diversity")


diss2 <-  ggplot(dd_dissimilarity_scenario2, aes(x = z1, y = rdiv, color = E1_mean)) +
  # scale_y_continuous(limits = c(0.99, 1.005), expand = c(0.02, 0.02)) +
  scale_color_uchicago() +
  labs(x = "Temperature optimum diversity", y = "Dissimilarity (derivatives)") +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 16),
    panel.grid = element_blank()
  ) +
  geom_jitter(size = 2, alpha = 0.25, width = 0.2) +
  stat_summary(fun = mean, geom = "point", size = 5) + labs(tag = "(b)")+ 
  scale_x_discrete(limits = level_order) + ggtitle("Treatment II:\nVariation in Temperature and Salinity\noptimum diversity - positive correlation")

diss3 <-  ggplot(dd_dissimilarity_scenario3, aes(x = z1, y = rdiv, color = E1_mean)) +
  # scale_y_continuous(limits = c(0.99, 1.005), expand = c(0.02, 0.02)) +
  scale_color_uchicago() +
  labs(x = "Temperature optimum diversity", y = "Dissimilarity (derivatives)") +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 16),
    panel.grid = element_blank()
  ) +
  geom_jitter(size = 2, alpha = 0.25, width = 0.2) +
  stat_summary(fun = mean, geom = "point", size = 5) + labs(tag = "(c)")+ 
  scale_x_discrete(limits = level_order) + ggtitle("Treatment II:\nVariation in Temperature and Salinity\noptimum diversity - negative correlation")


diss4 <-  ggplot(dd_dissimilarity_scenario1_SAME, aes(x = z1, y = rdiv, color = E1_mean)) +
  # scale_y_continuous(limits = c(0.99, 1.005), expand = c(0.02, 0.02)) +
  scale_color_uchicago() +
  labs(x = "Temperature optimum diversity", y = "Dissimilarity (derivatives)") +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 16),
    panel.grid = element_blank()
  ) +
  geom_jitter(size = 2, alpha = 0.25, width = 0.2) +
  stat_summary(fun = mean, geom = "point", size = 5) + labs(tag = "(d)")+ 
  scale_x_discrete(limits = level_order)  + ggtitle("Treatment II:\nOnly variation in Temperature\n optimum diversity")


diss5 <-  ggplot(dd_dissimilarity_scenario2_SAME, aes(x = z1, y = rdiv, color = E1_mean)) +
  # scale_y_continuous(limits = c(0.99, 1.005), expand = c(0.02, 0.02)) +
  scale_color_uchicago() +
  labs(x = "Temperature optimum diversity", y = "Dissimilarity (derivatives)") +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 16),
    panel.grid = element_blank()
  ) +
  geom_jitter(size = 2, alpha = 0.25, width = 0.2) +
  stat_summary(fun = mean, geom = "point", size = 5) + labs(tag = "(e)")+ 
  scale_x_discrete(limits = level_order) + ggtitle("Treatment II:\nVariation in Temperature and Salinity\noptimum diversity - positive correlation")

diss6 <-  ggplot(dd_dissimilarity_scenario3_SAME, aes(x = z1, y = rdiv, color = E1_mean)) +
  # scale_y_continuous(limits = c(0.99, 1.005), expand = c(0.02, 0.02)) +
  scale_color_uchicago() +
  labs(x = "Temperature optimum diversity", y = "Dissimilarity (derivatives)") +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 16),
    panel.grid = element_blank()
  ) +
  geom_jitter(size = 2, alpha = 0.25, width = 0.2) +
  stat_summary(fun = mean, geom = "point", size = 5) + labs(tag = "(f)")+ 
  scale_x_discrete(limits = level_order) + ggtitle("Treatment II:\nVariation in Temperature and Salinity\noptimum diversity - negative correlation")



top_plot      <- wrap_elements((diss1 + diss2 + diss3 & theme(legend.position = "bottom") & 
                                  labs(col = "Treatment III: Temperature mean value")) + plot_layout(guides = "collect")+
                                 plot_annotation(title = "Treatment I: Temperature dominant environmental variable",
                                                 theme = theme(plot.title = element_text(hjust = 0.5))))
  

bottom_plot   <- wrap_elements((diss4 + diss5 + diss6 & theme(legend.position = "bottom") & 
                                  labs(col = "Treatment III: Temperature mean value")) + plot_layout(guides = "collect")+ 
                                 plot_annotation(title = "Treatment I: Equal effect of Temperature and Salinity",
                                                 theme = theme(plot.title = element_text(hjust = 0.5))))

combined_final <- (top_plot / bottom_plot) 

combined_final <-combined_final + plot_layout(guides = "collect") # doesn't actually collect the legend... if you happen to know how to make it work, please get in touch
combined_final 

ggsave(combined_final, file = here("ms_figures/Fig12.jpg"), width = 18, height = 14)
ggsave(combined_final, file = here("ms_figures/Fig12.png"), width = 18, height = 14)
ggsave(combined_final, file = here("ms_figures/Fig12.pdf"), width = 18, height = 14)



############ Figure 13 ################
# Figure 13
theme_set(theme_classic(base_size = 15))
div1 <-  ggplot(dd_divergence_scenario1, aes(x = z1, y = sign, color = E1_mean)) +
  # scale_y_continuous(limits = c(0.99, 1.005), expand = c(0.02, 0.02)) +
  scale_color_uchicago() +
  labs(x = "Temperature optimum diversity", y = "Divergence (sign sensitive)") +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 16),
    panel.grid = element_blank()
  ) +
  geom_jitter(size = 2, alpha = 0.25, width = 0.2) +
  stat_summary(fun = mean, geom = "point", size = 5) + labs(tag = "(a)")+ 
  scale_x_discrete(limits = level_order)  + ggtitle("Treatment II:\nOnly variation in Temperature\n optimum diversity")


div2 <-  ggplot(dd_divergence_scenario2, aes(x = z1, y = sign, color = E1_mean)) +
  # scale_y_continuous(limits = c(0.99, 1.005), expand = c(0.02, 0.02)) +
  scale_color_uchicago() +
  labs(x = "Temperature optimum diversity", y = "Divergence (sign sensitive)") +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 16),
    panel.grid = element_blank()
  ) +
  geom_jitter(size = 2, alpha = 0.25, width = 0.2) +
  stat_summary(fun = mean, geom = "point", size = 5) + labs(tag = "(b)")+ 
  scale_x_discrete(limits = level_order) + ggtitle("Treatment II:\nVariation in Temperature and Salinity\noptimum diversity - positive correlation")

div3 <-  ggplot(dd_divergence_scenario3, aes(x = z1, y = sign, color = E1_mean)) +
  # scale_y_continuous(limits = c(0.99, 1.005), expand = c(0.02, 0.02)) +
  scale_color_uchicago() +
  labs(x = "Temperature optimum diversity", y = "Divergence (sign sensitive)") +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 16),
    panel.grid = element_blank()
  ) +
  geom_jitter(size = 2, alpha = 0.25, width = 0.2) +
  stat_summary(fun = mean, geom = "point", size = 5) + labs(tag = "(c)")+ 
  scale_x_discrete(limits = level_order) + ggtitle("Treatment II:\nVariation in Temperature and Salinity\noptimum diversity - negative correlation")


div4 <-  ggplot(dd_divergence_scenario1_SAME, aes(x = z1, y = sign, color = E1_mean)) +
  # scale_y_continuous(limits = c(0.99, 1.005), expand = c(0.02, 0.02)) +
  scale_color_uchicago() +
  labs(x = "Temperature optimum diversity", y = "Divergence (sign sensitive)") +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 16),
    panel.grid = element_blank()
  ) +
  geom_jitter(size = 2, alpha = 0.25, width = 0.2) +
  stat_summary(fun = mean, geom = "point", size = 5) + labs(tag = "(d)")+ 
  scale_x_discrete(limits = level_order)   + ggtitle("Treatment II:\nOnly variation in Temperature\n optimum diversity")


div5 <-  ggplot(dd_divergence_scenario2_SAME, aes(x = z1, y = sign, color = E1_mean)) +
  # scale_y_continuous(limits = c(0.99, 1.005), expand = c(0.02, 0.02)) +
  scale_color_uchicago() +
  labs(x = "Temperature optimum diversity", y = "Divergence (sign sensitive)") +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 16),
    panel.grid = element_blank()
  ) +
  geom_jitter(size = 2, alpha = 0.25, width = 0.2) +
  stat_summary(fun = mean, geom = "point", size = 5) + labs(tag = "(e)")+ 
  scale_x_discrete(limits = level_order) + ggtitle("Treatment II:\nVariation in Temperature and Salinity\noptimum diversity - positive correlation")


div6 <-  ggplot(dd_divergence_scenario3_SAME, aes(x = z1, y = sign, color = E1_mean)) +
  # scale_y_continuous(limits = c(0.99, 1.005), expand = c(0.02, 0.02)) +
  scale_color_uchicago() +
  labs(x = "Temperature optimum diversity", y = "Divergence (sign sensitive)") +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 16),
    panel.grid = element_blank()
  ) +
  geom_jitter(size = 2, alpha = 0.25, width = 0.2) +
  stat_summary(fun = mean, geom = "point", size = 5) + labs(tag = "(f)")+ 
  scale_x_discrete(limits = level_order) + ggtitle("Treatment II:\nVariation in Temperature and Salinity\noptimum diversity - negative correlation")





top_plot      <- wrap_elements((div1 + div2 + div3 & theme(legend.position = "bottom") & 
                                  labs(col = "Treatment III: Temperature mean value")) + plot_layout(guides = "collect")+
                                 plot_annotation(title = "Treatment I: Temperature dominant environmental variable",
                                                 theme = theme(plot.title = element_text(hjust = 0.5))))


bottom_plot   <- wrap_elements((div4 + div5 + div6 & theme(legend.position = "bottom") & 
                                  labs(col = "Treatment III: Temperature mean value")) + plot_layout(guides = "collect")+ 
                                 plot_annotation(title = "Treatment I: Equal effect of Temperature and Salinity",
                                                 theme = theme(plot.title = element_text(hjust = 0.5))))

combined_final <- (top_plot / bottom_plot) 

combined_final <-combined_final + plot_layout(guides = "collect") # doesn't actually collect the legend... if you happen to know how to meke it work, please get in touch
combined_final 

ggsave(combined_final, file = here("ms_figures/Fig13.jpg"), width = 18, height = 14)
ggsave(combined_final, file = here("ms_figures/Fig13.png"), width = 18, height = 14)
ggsave(combined_final, file = here("ms_figures/Fig13.pdf"), width = 18, height = 14)
