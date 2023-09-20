rm(list=ls())
library(tidyverse)
library(broom)
library(purrr)
library(mgcv)
library(gratia)
library(readr)
library(papeR)
library(here)
library(gtable)
library(gridExtra)
library(grid)
library(writexl)
library(plotly)
library(patchwork)
library(rlang)
library(plotly)
library(reshape2)
library(vctrs)
library(patchwork)
library(plot3D)
library(scales)
library(ggtext)
list.files("/Users/francesco/Documents/GitHub/multifarious_response_diversity/r/", full.names = TRUE) %>% map(source)
list.files("/Users/francesco/Documents/GitHub/gratia/R/", full.names = TRUE) %>% map(source)
list.files("/Users/francesco/Documents/GitHub/response_diversity_how_to_measure-main/r/", full.names = TRUE) %>% map(source)


## define the series of values of the environmental variables
E1_series <- seq(0, 50, 1)
E2_series <- seq(0, 50, 1)

## define new_data for following GAM fitting

new_data <- tibble(E1 = seq(0, 50, 0.2),
                   E2 = seq(0, 50, 0.2))





# Specify the correlation E1 and E2 
#negative in this example and high variation in E1 and E2 
theme_set(theme_classic(base_size = 12))
library(MASS)
sample_size <- 15                                       
sample_meanvector <- c(mean(E1_series),  mean(E2_series))                                   
sample_covariance_matrix <- matrix(c(150, -150,
                                     -150, 150),
                                   ncol = 2)

# create bivariate normal distribution
set.seed(65354)
refs <- mvrnorm(n = sample_size,
                mu = sample_meanvector, 
                Sigma = sample_covariance_matrix) %>%
  as_tibble() %>% 
  dplyr::rename(E1 = "V1", E2 = "V2") %>% 
  dplyr::mutate(time = seq.int(nrow(.)))

#cor(refs)

# plot env change over time
p_E1 <- refs %>% 
  ggplot(aes(x = time, y = E1)) + geom_line() +
  labs(tag = "(a)")

p_E2 <- refs %>% 
  ggplot(aes(x = time, y = E2)) + geom_line() +
  labs(tag = "(b)")

# Plot environmental change over time
p_E1 + p_E2
detach("package:MASS", unload=TRUE)

### Community 1 - low diversity in z1. 


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





s <- 4
# community 1 - low variation in z1. 
species_pars$z1_range <- 1.5
species_pars$z2_range <- 15

N_com <- seq(1, 10, 1)

# create reproducible community - 1
set.seed(1651)

# creating 20 communities with the same Z1_range and z2_range
communities1 <- rerun(10, get_community(s, species_pars))

# add community ID  communities 1 

mynames <-  as.list(N_com)

communities1 <-map2(communities1, mynames, ~.x %>% mutate(community = .y))

(communities1_tibble <-  tibble(
  sp = map(communities1, "species"),
  E1 = map(communities1, "E1"),
  E2 = map(communities1, "E2"),
  rate = map(communities1, "rate"),
  community = map(communities1, "community")) %>%
    unnest(sp, E1, E2, rate, community))


these_E2 <- seq( min(communities1_tibble$E1), max(communities1_tibble$E2), length.out= 6)
temp <- communities1_tibble %>%
  filter(E2 %in% these_E2) %>%
  mutate(E2 = as.factor(E2))
fig1 <- temp %>% 
  ggplot( aes(x = E1, y = rate, group = as.factor(E2), col = as.factor(E2))) +
  geom_point()+
  geom_smooth(se = FALSE) +
  facet_grid(community ~ sp) +
  theme_bw()+
  scale_color_viridis_d(end = 0.2,begin = 0.9, option = 'inferno') +
  labs(x = "E1", y = "rate", tag = "(a)")

these_E1 <- seq( min(communities1_tibble$E1), max(communities1_tibble$E1), length.out= 6)
temp <- communities1_tibble %>%
  filter(E1 %in% these_E1) %>%
  mutate(E1 = as.factor(E1))
fig2 <- temp %>% 
  ggplot( aes(x = E2, y = rate, group = as.factor(E1), col = as.factor(E1))) +
  geom_point()+
  geom_smooth(se = FALSE) +
  facet_grid(community ~ sp) +
  theme_bw()+
  scale_color_viridis_d(end = 0.2,begin = 0.9, option = 'inferno') +
  labs(x = "E2", y = "rate", tag = "(b)")

performance_comm1 <- fig1 + fig2
performance_comm1

# get directional derivatives 
dir_1 <- map(communities1, ~ .x %>% 
              get_directionals(., new_data, refs))


mynames <-  as.list(N_com)

dir_1 <-map2(dir_1, mynames, ~.x %>% mutate(community = .y))

(dir_1.2 <- tibble(
  sp = map(dir_1, "sp"),
  time = map(dir_1, "time"),
  E1 = map(dir_1, "E1_ref"),
  E2 = map(dir_1, "E2_ref"),
  dir_deriv = map(dir_1, "dir_deriv"),
  community = map(dir_1, "community")) %>%
    unnest(sp, time, E1, E1, dir_deriv, community))

# from long to wide
rdiv_1 <- map(dir_1, ~ .x %>% 
             spread(key = sp, value = dir_deriv))


### Response diversity calculation
# dissimilarity

rdiv_1 <- lapply(rdiv_1, function(x) cbind(x,"rdiv"= apply(subset(x, select = c(s1, s2, s3, s4)), 1, resp_div, sign_sens = F)))
rdiv_1 <- lapply(rdiv_1, function(x) cbind(x,"Med"= apply(subset(x, select = c(rdiv)),2, median)))

# divergence
rdiv_1 <- lapply(rdiv_1, function(x) cbind(x,"sign"= apply(subset(x, select = c(s1, s2, s3, s4)), 1, resp_div, sign_sens = T)))
rdiv_1 <- lapply(rdiv_1, function(x) cbind(x,"Med_sing"= apply(subset(x, select = c(sign)),2, median)))


# from list to tibble
(rdiv_1.2 <- tibble(
  time = map(rdiv_1, "time"),
  E1 = map(rdiv_1, "E1_ref"),
  E2 = map(rdiv_1, "E2_ref"),
  rdiv = map(rdiv_1, "rdiv"),
  sign = map(rdiv_1, "sign"),
  community = map(rdiv_1, "community")) %>%
    unnest(time, E1, E2, rdiv, sign, community))

rdiv_1.2 <- rdiv_1.2 %>%  group_by(community) %>% 
  mutate(Med = mean(rdiv), 
         Med_sign = mean(sign))

# plotting 
rmax<-max(rdiv_1.2$rdiv); rmin<-1
dmax<-max(rdiv_1.2$sign); dmin<-0

Fig_1_comm1 <-dir_1.2 %>%
  ggplot(aes(x=time, y=dir_deriv, col = sp)) +
  #theme_classic(base_size = 14) + 
  labs(x = "time",y = "Directional derivative",tag = "a)") + 
  facet_wrap(vars(community), nrow = 2) +
  #geom_point(size=0.5) +
  geom_line() + theme_bw() + geom_hline(yintercept = 0, linetype= "dashed")


Fig2_comm1 <- rdiv_1.2 %>% 
  ggplot(mapping = aes(x = time,y = rdiv)) +
  theme_bw(base_size = 12)+
  labs(x = "time",y = "Dissimilarity", tag = "b)") + 
  facet_wrap(vars(community), nrow = 2) +
  geom_line(aes(x=time, y=Med), colour="grey50", lty=2)+
  geom_line() + 
  geom_richtext(x = 7,
                mapping = aes(y = rmax),
                size=4.5,
                label.color = NA,
                label = paste("RDiv^Mean =",paste0(round(rdiv_1.2$Med,digits = 4))))+
  lims(y = c(rmin,rmax))



Fig3_comm1 <- rdiv_1.2 %>% 
  ggplot(mapping = aes(x = time,y = sign)) +
  theme_bw(base_size = 12)+
  labs(x = "time",y = "Divergence", tag = "c)") + 
  facet_wrap(vars(community), nrow = 2) +
  geom_line(aes(x=time, y=Med_sign), colour="grey50", lty=2)+
  geom_line() + 
  geom_richtext(x = 7,
                mapping = aes(y = dmax),
                size=4.5,
                label.color = NA,
                label = paste("RDiv^Mean =",paste0(round(rdiv_1.2$Med_sign,digits = 4))))+
  lims(y = c(dmin,dmax))

Fig_1_comm1 / Fig2_comm1/ Fig3_comm1


