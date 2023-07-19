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
E1_series <- 273.15 + seq(0, 50, 1)
E2_series <- seq(0, 50, 1)

## define new_data for following GAM fitting

new_data <- tibble(E1 = 273.15 + seq(0, 50, 0.2),
                   E2 = seq(0, 50, 0.2))


#######################################################
###### Variation in z1 only #########
#######################################################

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
species_pars <- list(a1_mean = 1e-9,
                     b1_mean = 0.063,
                     z1_mean = 285,
                     w1_mean = 60,
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


these_E2 <- seq( min(communities1_tibble$E2), max(communities1_tibble$E2), length.out= 6)
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




#### Medium variation in z1 - negative correlation E1 E2 
s <- 4
# community 1 - low variation in z1. 
species_pars$z1_range <- 25
species_pars$z2_range <- 15

N_com <- seq(1, 10, 1)

# create reproducible community - 1
set.seed(1651)

# creating 20 communities with the same Z1_range and z2_range
communities2 <- rerun(10, get_community(s, species_pars))

# add community ID  communities 1 

mynames <-  as.list(N_com)

communities2 <-map2(communities2, mynames, ~.x %>% mutate(community = .y))

(communities2_tibble <-  tibble(
  sp = map(communities2, "species"),
  E1 = map(communities2, "E1"),
  E2 = map(communities2, "E2"),
  rate = map(communities2, "rate"),
  community = map(communities2, "community")) %>%
    unnest(sp, E1, E2, rate, community))


these_E2 <- seq( min(communities2_tibble$E2), max(communities2_tibble$E2), length.out= 6)
temp <- communities2_tibble %>%
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

these_E1 <- seq( min(communities2_tibble$E1), max(communities2_tibble$E1), length.out= 6)
temp <- communities2_tibble %>%
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

performance_comm2 <- fig1 + fig2
performance_comm2

# get directional derivatives 
dir_2 <- map(communities2, ~ .x %>% 
               get_directionals(., new_data, refs))


mynames <-  as.list(N_com)

dir_2 <-map2(dir_2, mynames, ~.x %>% mutate(community = .y))

(dir_2.2 <- tibble(
  sp = map(dir_2, "sp"),
  time = map(dir_2, "time"),
  E1 = map(dir_2, "E1_ref"),
  E2 = map(dir_2, "E2_ref"),
  dir_deriv = map(dir_2, "dir_deriv"),
  community = map(dir_2, "community")) %>%
    unnest(sp, time, E1, E1, dir_deriv, community))

# from long to wide
rdiv_2 <- map(dir_2, ~ .x %>% 
                spread(key = sp, value = dir_deriv))


### Response diversity calculation
# dissimilarity
rdiv_2 <- lapply(rdiv_2, function(x) cbind(x,"rdiv"= apply(subset(x, select = c(s1, s2, s3, s4)), 1, resp_div, sign_sens = F)))
rdiv_2 <- lapply(rdiv_2, function(x) cbind(x,"Med"= apply(subset(x, select = c(rdiv)),2, median)))

# divergence
rdiv_2 <- lapply(rdiv_2, function(x) cbind(x,"sign"= apply(subset(x, select = c(s1, s2, s3, s4)), 1, resp_div, sign_sens = T)))
rdiv_2 <- lapply(rdiv_2, function(x) cbind(x,"Med_sing"= apply(subset(x, select = c(sign)),2, median)))


# from list to tibble
(rdiv_2.2 <- tibble(
  time = map(rdiv_2, "time"),
  E1 = map(rdiv_2, "E1_ref"),
  E2 = map(rdiv_2, "E2_ref"),
  rdiv = map(rdiv_2, "rdiv"),
  sign = map(rdiv_2, "sign"),
  community = map(rdiv_2, "community")) %>%
    unnest(time, E1, E2, rdiv, sign, community))

rdiv_2.2 <- rdiv_2.2 %>%  group_by(community) %>% 
  mutate(Med = mean(rdiv), 
         Med_sign = mean(sign))

# plotting 
rmax<-max(rdiv_2.2$rdiv); rmin<-1
dmax<-max(rdiv_2.2$sign); dmin<-0

Fig_1_comm2 <-dir_2.2 %>%
  ggplot(aes(x=time, y=dir_deriv, col = sp)) +
  #theme_classic(base_size = 14) + 
  labs(x = "time",y = "Directional derivative",tag = "a)") + 
  facet_wrap(vars(community), nrow = 2) +
  #geom_point(size=0.5) +
  geom_line() + theme_bw() + geom_hline(yintercept = 0, linetype= "dashed")


Fig2_comm2 <- rdiv_2.2 %>% 
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
                label = paste("RDiv^Mean =",paste0(round(rdiv_2.2$Med,digits = 4))))+
  lims(y = c(rmin,rmax))



Fig3_comm2 <- rdiv_2.2 %>% 
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
                label = paste("RDiv^Mean =",paste0(round(rdiv_2.2$Med_sign,digits = 4))))+
  lims(y = c(dmin,dmax))

Fig_1_comm2 / Fig2_comm2/ Fig3_comm2





#### High variation in z1 - negative correlation E1 E2 
s <- 4
# community 1 - low variation in z1. 
species_pars$z1_range <- 50
species_pars$z2_range <- 15

N_com <- seq(1, 10, 1)

# create reproducible community - 1
set.seed(1651)

# creating 20 communities with the same Z1_range and z2_range
communities3 <- rerun(10, get_community(s, species_pars))

# add community ID  communities 1 

mynames <-  as.list(N_com)

communities3 <-map2(communities3, mynames, ~.x %>% mutate(community = .y))

(communities3_tibble <-  tibble(
  sp = map(communities3, "species"),
  E1 = map(communities3, "E1"),
  E2 = map(communities3, "E2"),
  rate = map(communities3, "rate"),
  community = map(communities3, "community")) %>%
    unnest(sp, E1, E2, rate, community))


these_E2 <- seq( min(communities3_tibble$E2), max(communities3_tibble$E2), length.out= 6)
temp <- communities3_tibble %>%
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

these_E1 <- seq( min(communities3_tibble$E1), max(communities3_tibble$E1), length.out= 6)
temp <- communities3_tibble %>%
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

performance_comm3 <- fig1 + fig2
performance_comm3

# get directional derivatives 
dir_3 <- map(communities3, ~ .x %>% 
               get_directionals(., new_data, refs))


mynames <-  as.list(N_com)

dir_3 <-map2(dir_3, mynames, ~.x %>% mutate(community = .y))

(dir_3.2 <- tibble(
  sp = map(dir_3, "sp"),
  time = map(dir_3, "time"),
  E1 = map(dir_3, "E1_ref"),
  E2 = map(dir_3, "E2_ref"),
  dir_deriv = map(dir_3, "dir_deriv"),
  community = map(dir_3, "community")) %>%
    unnest(sp, time, E1, E1, dir_deriv, community))

# from long to wide
rdiv_3 <- map(dir_3, ~ .x %>% 
                spread(key = sp, value = dir_deriv))


### Response diversity calculation
# dissimilarity
rdiv_3 <- lapply(rdiv_3, function(x) cbind(x,"rdiv"= apply(subset(x, select = c(s1, s2, s3, s4)), 1, resp_div, sign_sens = F)))
rdiv_3 <- lapply(rdiv_3, function(x) cbind(x,"Med"= apply(subset(x, select = c(rdiv)),2, median)))

# divergence
rdiv_3 <- lapply(rdiv_3, function(x) cbind(x,"sign"= apply(subset(x, select = c(s1, s2, s3, s4)), 1, resp_div, sign_sens = T)))
rdiv_3 <- lapply(rdiv_3, function(x) cbind(x,"Med_sing"= apply(subset(x, select = c(sign)),2, median)))


# from list to tibble
(rdiv_3.2 <- tibble(
  time = map(rdiv_3, "time"),
  E1 = map(rdiv_3, "E1_ref"),
  E2 = map(rdiv_3, "E2_ref"),
  rdiv = map(rdiv_3, "rdiv"),
  sign = map(rdiv_3, "sign"),
  community = map(rdiv_3, "community")) %>%
    unnest(time, E1, E2, rdiv, sign, community))

rdiv_3.2 <- rdiv_3.2 %>%  group_by(community) %>% 
  mutate(Med = mean(rdiv), 
         Med_sign = mean(sign))

# plotting 
rmax<-max(rdiv_3.2$rdiv); rmin<-1
dmax<-max(rdiv_3.2$sign); dmin<-0

Fig_1_comm3 <-dir_3.2 %>%
  ggplot(aes(x=time, y=dir_deriv, col = sp)) +
  #theme_classic(base_size = 14) + 
  labs(x = "time",y = "Directional derivative",tag = "a)") + 
  facet_wrap(vars(community), nrow = 2) +
  #geom_point(size=0.5) +
  geom_line() + theme_bw() + geom_hline(yintercept = 0, linetype= "dashed")


Fig2_comm3 <- rdiv_3.2 %>% 
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
                label = paste("RDiv^Mean =",paste0(round(rdiv_2.2$Med,digits = 4))))+
  lims(y = c(rmin,rmax))



Fig3_comm3 <- rdiv_3.2 %>% 
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
                label = paste("RDiv^Mean =",paste0(round(rdiv_2.2$Med_sign,digits = 4))))+
  lims(y = c(dmin,dmax))

Fig_1_comm3 / Fig2_comm3/ Fig3_comm3

### For summary plot

dd1 <- data.frame(rdiv_1.2[, c(4:6)]) 
dd1$cor <- "negative"
dd1$z1 <- "low"

dd2 <- data.frame(rdiv_2.2[, c(4:6)]) 
dd2$cor <- "negative"
dd2$z1 <- "medium"

dd3 <- data.frame(rdiv_3.2[, c(4:6)]) 
dd3$cor <- "negative"
dd3$z1 <- "high"




## High fluctuation in environmental change - positive correlation
#Increasing variance in z1 - positive correlation between env variables (E1 and E2) with *high* fluctuations

#Same steps as before (3 communities with increasing diversity ib z1 while z2 is fixed), but the correlation between E1 and E1 is negative.

library(MASS)
sample_size <- 15                                       
sample_medianvector <- c(mean(E1_series),  mean(E2_series))                                   
sample_covariance_matrix <- matrix(c(150, 150,
                                     150, 150),
                                   ncol = 2)

# create bivariate normal distribution
set.seed(5655)
refs <- mvrnorm(n = sample_size,
                mu = sample_meanvector, 
                Sigma = sample_covariance_matrix) %>%
  as_tibble() %>% 
  dplyr::rename(E1 = "V1", E2 = "V2") %>% 
  dplyr::mutate(time = seq.int(nrow(.)))



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



### Creating communities

# community 1 
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
# get directional derivatives 
dir_1 <- map(communities1, ~ .x %>% 
               get_directionals(., new_data, refs))

mynames <-  as.list(N_com)

dir_1 <-map2(dir_1, mynames, ~.x %>% mutate(community = .y))


# community 2 
s <- 4
# community 1 - low variation in z1. 
species_pars$z1_range <- 25
species_pars$z2_range <- 15

N_com <- seq(1, 10, 1)

# create reproducible community - 2
set.seed(1651)

# creating 20 communities with the same Z1_range and z2_range
communities2 <- rerun(10, get_community(s, species_pars))

# add community ID  communities 2 

mynames <-  as.list(N_com)

communities2 <-map2(communities2, mynames, ~.x %>% mutate(community = .y))
# get directional derivatives 
dir_2 <- map(communities2, ~ .x %>% 
               get_directionals(., new_data, refs))


mynames <-  as.list(N_com)

dir_2 <-map2(dir_2, mynames, ~.x %>% mutate(community = .y))

# community 3
s <- 4
# community 1 - low variation in z1. 
species_pars$z1_range <- 50
species_pars$z2_range <- 15

N_com <- seq(1, 10, 1)

# create reproducible community - 1
set.seed(1651)

# creating 20 communities with the same Z1_range and z2_range
communities3 <- rerun(10, get_community(s, species_pars))

# add community ID  communities 1 

mynames <-  as.list(N_com)

communities3 <-map2(communities3, mynames, ~.x %>% mutate(community = .y))
# get directional derivatives 
dir_3 <- map(communities3, ~ .x %>% 
               get_directionals(., new_data, refs))


mynames <-  as.list(N_com)

dir_3 <-map2(dir_3, mynames, ~.x %>% mutate(community = .y))


### Response diversity calculation

#comm 1 
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



Fig3_comm1<- rdiv_1.2 %>% 
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


#comm 2
# from long to wide
rdiv_2 <- map(dir_2, ~ .x %>% 
                spread(key = sp, value = dir_deriv))


### Response diversity calculation
# dissimilarity
rdiv_2 <- lapply(rdiv_2, function(x) cbind(x,"rdiv"= apply(subset(x, select = c(s1, s2, s3, s4)), 1, resp_div, sign_sens = F)))
rdiv_2 <- lapply(rdiv_2, function(x) cbind(x,"Med"= apply(subset(x, select = c(rdiv)),2, median)))

# divergence
rdiv_2 <- lapply(rdiv_2, function(x) cbind(x,"sign"= apply(subset(x, select = c(s1, s2, s3, s4)), 1, resp_div, sign_sens = T)))
rdiv_2 <- lapply(rdiv_2, function(x) cbind(x,"Med_sing"= apply(subset(x, select = c(sign)),2, median)))


# from list to tibble
(rdiv_2.2 <- tibble(
  time = map(rdiv_2, "time"),
  E1 = map(rdiv_2, "E1_ref"),
  E2 = map(rdiv_2, "E2_ref"),
  rdiv = map(rdiv_2, "rdiv"),
  sign = map(rdiv_2, "sign"),
  community = map(rdiv_2, "community")) %>%
    unnest(time, E1, E2, rdiv, sign, community))

rdiv_2.2 <- rdiv_2.2 %>%  group_by(community) %>% 
  mutate(Med = mean(rdiv), 
         Med_sign = mean(sign))

# plotting 
rmax<-max(rdiv_2.2$rdiv); rmin<-1
dmax<-max(rdiv_2.2$sign); dmin<-0

Fig_1_comm2 <-dir_2.2 %>%
  ggplot(aes(x=time, y=dir_deriv, col = sp)) +
  #theme_classic(base_size = 14) + 
  labs(x = "time",y = "Directional derivative",tag = "a)") + 
  facet_wrap(vars(community), nrow = 2) +
  #geom_point(size=0.5) +
  geom_line() + theme_bw() + geom_hline(yintercept = 0, linetype= "dashed")


Fig2_comm2 <- rdiv_2.2 %>% 
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
                label = paste("RDiv^Mean =",paste0(round(rdiv_2.2$Med,digits = 4))))+
  lims(y = c(rmin,rmax))



Fig3_comm2<- rdiv_2.2 %>% 
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
                label = paste("RDiv^Mean =",paste0(round(rdiv_2.2$Med_sign,digits = 4))))+
  lims(y = c(dmin,dmax))

Fig_1_comm2 / Fig2_comm2/ Fig3_comm2




#comm 3
# from long to wide
rdiv_3 <- map(dir_3, ~ .x %>% 
                spread(key = sp, value = dir_deriv))


### Response diversity calculation
# dissimilarity
rdiv_3 <- lapply(rdiv_3, function(x) cbind(x,"rdiv"= apply(subset(x, select = c(s1, s2, s3, s4)), 1, resp_div, sign_sens = F)))
rdiv_3 <- lapply(rdiv_3, function(x) cbind(x,"Med"= apply(subset(x, select = c(rdiv)),2, median)))

# divergence
rdiv_3 <- lapply(rdiv_3, function(x) cbind(x,"sign"= apply(subset(x, select = c(s1, s2, s3, s4)), 1, resp_div, sign_sens = T)))
rdiv_3 <- lapply(rdiv_3, function(x) cbind(x,"Med_sing"= apply(subset(x, select = c(sign)),2, median)))


# from list to tibble
(rdiv_3.2 <- tibble(
  time = map(rdiv_3, "time"),
  E1 = map(rdiv_3, "E1_ref"),
  E2 = map(rdiv_3, "E2_ref"),
  rdiv = map(rdiv_3, "rdiv"),
  sign = map(rdiv_3, "sign"),
  community = map(rdiv_3, "community")) %>%
    unnest(time, E1, E2, rdiv, sign, community))

rdiv_3.2 <- rdiv_3.2 %>%  group_by(community) %>% 
  mutate(Med = mean(rdiv), 
         Med_sign = mean(sign))

# plotting 
rmax<-max(rdiv_3.2$rdiv); rmin<-1
dmax<-max(rdiv_3.2$sign); dmin<-0

Fig_1_comm3 <-dir_3.2 %>%
  ggplot(aes(x=time, y=dir_deriv, col = sp)) +
  #theme_classic(base_size = 14) + 
  labs(x = "time",y = "Directional derivative",tag = "a)") + 
  facet_wrap(vars(community), nrow = 2) +
  #geom_point(size=0.5) +
  geom_line() + theme_bw() + geom_hline(yintercept = 0, linetype= "dashed")


Fig2_comm3 <- rdiv_3.2 %>% 
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
                label = paste("RDiv^Mean =",paste0(round(rdiv_3.2$Med,digits = 4))))+
  lims(y = c(rmin,rmax))



Fig3_comm3<- rdiv_3.2 %>% 
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
                label = paste("RDiv^Mean =",paste0(round(rdiv_3.2$Med_sign,digits = 4))))+
  lims(y = c(dmin,dmax))

Fig_1_comm3/ Fig2_comm3/ Fig3_comm3


### For summary plot

dd4 <- data.frame(rdiv_1.2[, c(4:6)]) 
dd4$cor <- "positive"
dd4$z1 <- "low"

dd5 <- data.frame(rdiv_2.2[, c(4:6)]) 
dd5$cor <- "positive"
dd5$z1 <- "medium"

dd6 <- data.frame(rdiv_3.2[, c(4:6)]) 
dd6$cor <- "positive"
dd6$z1 <- "high"




## High fluctuation in environmental change - no correlation
# Increasing variance in z1 - no correlation between env variables (E1 and E2) with *high* fluctuations
# 
# Same steps as before (3 communities with increasing diversity in z1 while z2 is fixed), but there is no correlation between E1 and E1

library(MASS)
sample_size <- 15                                       
sample_meanvector <- c(mean(E1_series),  mean(E2_series))                                   
sample_covariance_matrix <- matrix(c(150, 0.0,
                                     0.0, 90),
                                   ncol = 2)

# create bivariate normal distribution
set.seed(5655)
refs <- mvrnorm(n = sample_size,
                mu = sample_meanvector, 
                Sigma = sample_covariance_matrix) %>%
  as_tibble() %>% 
  dplyr::rename(E1 = "V1", E2 = "V2") %>% 
  dplyr::mutate(time = seq.int(nrow(.)))


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






### Creating communities

# community 1 
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
# get directional derivatives 
dir_1 <- map(communities1, ~ .x %>% 
               get_directionals(., new_data, refs))

mynames <-  as.list(N_com)

dir_1 <-map2(dir_1, mynames, ~.x %>% mutate(community = .y))


# community 2 
s <- 4
# community 1 - low variation in z1. 
species_pars$z1_range <- 25
species_pars$z2_range <- 15

N_com <- seq(1, 10, 1)

# create reproducible community - 2
set.seed(1651)

# creating 20 communities with the same Z1_range and z2_range
communities2 <- rerun(10, get_community(s, species_pars))

# add community ID  communities 2 

mynames <-  as.list(N_com)

communities2 <-map2(communities2, mynames, ~.x %>% mutate(community = .y))
# get directional derivatives 
dir_2 <- map(communities2, ~ .x %>% 
               get_directionals(., new_data, refs))


mynames <-  as.list(N_com)

dir_2 <-map2(dir_2, mynames, ~.x %>% mutate(community = .y))

# community 3
s <- 4
# community 1 - low variation in z1. 
species_pars$z1_range <- 50
species_pars$z2_range <- 15

N_com <- seq(1, 10, 1)

# create reproducible community - 1
set.seed(1651)

# creating 20 communities with the same Z1_range and z2_range
communities3 <- rerun(10, get_community(s, species_pars))

# add community ID  communities 1 

mynames <-  as.list(N_com)

communities3 <-map2(communities3, mynames, ~.x %>% mutate(community = .y))
# get directional derivatives 
dir_3 <- map(communities3, ~ .x %>% 
               get_directionals(., new_data, refs))


mynames <-  as.list(N_com)

dir_3 <-map2(dir_3, mynames, ~.x %>% mutate(community = .y))


### Response diversity calculation

#comm 1 
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



Fig3_comm1<- rdiv_1.2 %>% 
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


#comm 2
# from long to wide
rdiv_2 <- map(dir_2, ~ .x %>% 
                spread(key = sp, value = dir_deriv))


### Response diversity calculation
# dissimilarity
rdiv_2 <- lapply(rdiv_2, function(x) cbind(x,"rdiv"= apply(subset(x, select = c(s1, s2, s3, s4)), 1, resp_div, sign_sens = F)))
rdiv_2 <- lapply(rdiv_2, function(x) cbind(x,"Med"= apply(subset(x, select = c(rdiv)),2, median)))

# divergence
rdiv_2 <- lapply(rdiv_2, function(x) cbind(x,"sign"= apply(subset(x, select = c(s1, s2, s3, s4)), 1, resp_div, sign_sens = T)))
rdiv_2 <- lapply(rdiv_2, function(x) cbind(x,"Med_sing"= apply(subset(x, select = c(sign)),2, median)))


# from list to tibble
(rdiv_2.2 <- tibble(
  time = map(rdiv_2, "time"),
  E1 = map(rdiv_2, "E1_ref"),
  E2 = map(rdiv_2, "E2_ref"),
  rdiv = map(rdiv_2, "rdiv"),
  sign = map(rdiv_2, "sign"),
  community = map(rdiv_2, "community")) %>%
    unnest(time, E1, E2, rdiv, sign, community))

rdiv_2.2 <- rdiv_2.2 %>%  group_by(community) %>% 
  mutate(Med = mean(rdiv), 
         Med_sign = mean(sign))

# plotting 
rmax<-max(rdiv_2.2$rdiv); rmin<-1
dmax<-max(rdiv_2.2$sign); dmin<-0

Fig_1_comm2 <-dir_2.2 %>%
  ggplot(aes(x=time, y=dir_deriv, col = sp)) +
  #theme_classic(base_size = 14) + 
  labs(x = "time",y = "Directional derivative",tag = "a)") + 
  facet_wrap(vars(community), nrow = 2) +
  #geom_point(size=0.5) +
  geom_line() + theme_bw() + geom_hline(yintercept = 0, linetype= "dashed")


Fig2_comm2 <- rdiv_2.2 %>% 
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
                label = paste("RDiv^Mean =",paste0(round(rdiv_2.2$Med,digits = 4))))+
  lims(y = c(rmin,rmax))



Fig3_comm2<- rdiv_2.2 %>% 
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
                label = paste("RDiv^Mean =",paste0(round(rdiv_2.2$Med_sign,digits = 4))))+
  lims(y = c(dmin,dmax))

Fig_1_comm2 / Fig2_comm2/ Fig3_comm2




#comm 3
# from long to wide
rdiv_3 <- map(dir_3, ~ .x %>% 
                spread(key = sp, value = dir_deriv))


### Response diversity calculation
# dissimilarity
rdiv_3 <- lapply(rdiv_3, function(x) cbind(x,"rdiv"= apply(subset(x, select = c(s1, s2, s3, s4)), 1, resp_div, sign_sens = F)))
rdiv_3 <- lapply(rdiv_3, function(x) cbind(x,"Med"= apply(subset(x, select = c(rdiv)),2, median)))

# divergence
rdiv_3 <- lapply(rdiv_3, function(x) cbind(x,"sign"= apply(subset(x, select = c(s1, s2, s3, s4)), 1, resp_div, sign_sens = T)))
rdiv_3 <- lapply(rdiv_3, function(x) cbind(x,"Med_sing"= apply(subset(x, select = c(sign)),2, median)))


# from list to tibble
(rdiv_3.2 <- tibble(
  time = map(rdiv_3, "time"),
  E1 = map(rdiv_3, "E1_ref"),
  E2 = map(rdiv_3, "E2_ref"),
  rdiv = map(rdiv_3, "rdiv"),
  sign = map(rdiv_3, "sign"),
  community = map(rdiv_3, "community")) %>%
    unnest(time, E1, E2, rdiv, sign, community))

rdiv_3.2 <- rdiv_3.2 %>%  group_by(community) %>% 
  mutate(Med = mean(rdiv), 
         Med_sign = mean(sign))

# plotting 
rmax<-max(rdiv_3.2$rdiv); rmin<-1
dmax<-max(rdiv_3.2$sign); dmin<-0

Fig_1_comm3 <-dir_3.2 %>%
  ggplot(aes(x=time, y=dir_deriv, col = sp)) +
  #theme_classic(base_size = 14) + 
  labs(x = "time",y = "Directional derivative",tag = "a)") + 
  facet_wrap(vars(community), nrow = 2) +
  #geom_point(size=0.5) +
  geom_line() + theme_bw() + geom_hline(yintercept = 0, linetype= "dashed")


Fig2_comm3 <- rdiv_3.2 %>% 
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
                label = paste("RDiv^Mean =",paste0(round(rdiv_3.2$Med,digits = 4))))+
  lims(y = c(rmin,rmax))



Fig3_comm3<- rdiv_3.2 %>% 
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
                label = paste("RDiv^Mean =",paste0(round(rdiv_3.2$Med_sign,digits = 4))))+
  lims(y = c(dmin,dmax))

Fig_1_comm3/ Fig2_comm3/ Fig3_comm3


### For summary plot

dd7 <- data.frame(rdiv_1.2[, c(4:6)]) 
dd7$cor <- "no_cor"
dd7$z1 <- "low"

dd8 <- data.frame(rdiv_2.2[, c(4:6)]) 
dd8$cor <- "no_cor"
dd8$z1 <- "medium"

dd9 <- data.frame(rdiv_3.2[, c(4:6)]) 
dd9$cor <- "no_cor"
dd9$z1 <- "high"



dd_plot1 <- rbind(dd1, dd2, dd3,
                  dd4, dd5, dd6, 
                  dd7, dd8, dd9)

theme_set(theme_classic(base_size = 12))

g1.1 <-
  ggplot(dd_plot1, aes(x = cor, y = rdiv, color = z1)) +
  # scale_y_continuous(limits = c(0.99, 1.005), expand = c(0.02, 0.02)) +
  scale_color_uchicago() +
  labs(x = NULL, y = "Dissimilarity (derivatives)") +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 16),
    panel.grid = element_blank()
  )


rdiv1 <- g1.1 +
  geom_jitter(size = 2, alpha = 0.25, width = 0.2) +
  stat_summary(fun = mean, geom = "point", size = 5) + labs(tag = "(a)")



g1.2 <-
  ggplot(dd_plot1, aes(x = cor, y = sign, color = z1)) +
  scale_color_uchicago() +
  labs(x = NULL, y = "Divergence (sign sensitive)") +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 16),
    panel.grid = element_blank()
  )

sign1 <- g1.2 +
  geom_jitter(size = 2, alpha = 0.25, width = 0.2) +
  stat_summary(fun = mean, geom = "point", size = 5)+ labs(tag = "(b)") +
  ylim(0, 1)

combined <- rdiv1 + sign1 & theme(legend.position = "bottom")
combined + plot_layout(guides = "collect")











#######################################################
###### Variation in z1 and z2 - both increase #########
#######################################################
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
species_pars <- list(a1_mean = 1e-9,
                     b1_mean = 0.063,
                     z1_mean = 285,
                     w1_mean = 60,
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
species_pars$z2_range <- 1.5

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


these_E2 <- seq( min(communities1_tibble$E2), max(communities1_tibble$E2), length.out= 6)
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




#### Medium variation in z1 and z2 - negative correlation E1 E2 
s <- 4
# community 1 - low variation in z1. 
species_pars$z1_range <- 25
species_pars$z2_range <- 25

N_com <- seq(1, 10, 1)

# create reproducible community - 1
set.seed(1651)

# creating 20 communities with the same Z1_range and z2_range
communities2 <- rerun(10, get_community(s, species_pars))

# add community ID  communities 1 

mynames <-  as.list(N_com)

communities2 <-map2(communities2, mynames, ~.x %>% mutate(community = .y))

(communities2_tibble <-  tibble(
  sp = map(communities2, "species"),
  E1 = map(communities2, "E1"),
  E2 = map(communities2, "E2"),
  rate = map(communities2, "rate"),
  community = map(communities2, "community")) %>%
    unnest(sp, E1, E2, rate, community))


these_E2 <- seq( min(communities2_tibble$E2), max(communities2_tibble$E2), length.out= 6)
temp <- communities2_tibble %>%
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

these_E1 <- seq( min(communities2_tibble$E1), max(communities2_tibble$E1), length.out= 6)
temp <- communities2_tibble %>%
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

performance_comm2 <- fig1 + fig2
performance_comm2

# get directional derivatives 
dir_2 <- map(communities2, ~ .x %>% 
               get_directionals(., new_data, refs))


mynames <-  as.list(N_com)

dir_2 <-map2(dir_2, mynames, ~.x %>% mutate(community = .y))

(dir_2.2 <- tibble(
  sp = map(dir_2, "sp"),
  time = map(dir_2, "time"),
  E1 = map(dir_2, "E1_ref"),
  E2 = map(dir_2, "E2_ref"),
  dir_deriv = map(dir_2, "dir_deriv"),
  community = map(dir_2, "community")) %>%
    unnest(sp, time, E1, E1, dir_deriv, community))

# from long to wide
rdiv_2 <- map(dir_2, ~ .x %>% 
                spread(key = sp, value = dir_deriv))


### Response diversity calculation
# dissimilarity
rdiv_2 <- lapply(rdiv_2, function(x) cbind(x,"rdiv"= apply(subset(x, select = c(s1, s2, s3, s4)), 1, resp_div, sign_sens = F)))
rdiv_2 <- lapply(rdiv_2, function(x) cbind(x,"Med"= apply(subset(x, select = c(rdiv)),2, median)))

# divergence
rdiv_2 <- lapply(rdiv_2, function(x) cbind(x,"sign"= apply(subset(x, select = c(s1, s2, s3, s4)), 1, resp_div, sign_sens = T)))
rdiv_2 <- lapply(rdiv_2, function(x) cbind(x,"Med_sing"= apply(subset(x, select = c(sign)),2, median)))


# from list to tibble
(rdiv_2.2 <- tibble(
  time = map(rdiv_2, "time"),
  E1 = map(rdiv_2, "E1_ref"),
  E2 = map(rdiv_2, "E2_ref"),
  rdiv = map(rdiv_2, "rdiv"),
  sign = map(rdiv_2, "sign"),
  community = map(rdiv_2, "community")) %>%
    unnest(time, E1, E2, rdiv, sign, community))

rdiv_2.2 <- rdiv_2.2 %>%  group_by(community) %>% 
  mutate(Med = mean(rdiv), 
         Med_sign = mean(sign))

# plotting 
rmax<-max(rdiv_2.2$rdiv); rmin<-1
dmax<-max(rdiv_2.2$sign); dmin<-0

Fig_1_comm2 <-dir_2.2 %>%
  ggplot(aes(x=time, y=dir_deriv, col = sp)) +
  #theme_classic(base_size = 14) + 
  labs(x = "time",y = "Directional derivative",tag = "a)") + 
  facet_wrap(vars(community), nrow = 2) +
  #geom_point(size=0.5) +
  geom_line() + theme_bw() + geom_hline(yintercept = 0, linetype= "dashed")


Fig2_comm2 <- rdiv_2.2 %>% 
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
                label = paste("RDiv^Mean =",paste0(round(rdiv_2.2$Med,digits = 4))))+
  lims(y = c(rmin,rmax))



Fig3_comm2 <- rdiv_2.2 %>% 
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
                label = paste("RDiv^Mean =",paste0(round(rdiv_2.2$Med_sign,digits = 4))))+
  lims(y = c(dmin,dmax))

Fig_1_comm2 / Fig2_comm2/ Fig3_comm2





#### High variation in z1 - negative correlation E1 E2 
s <- 4
# community 1 - low variation in z1. 
species_pars$z1_range <- 50
species_pars$z2_range <- 50

N_com <- seq(1, 10, 1)

# create reproducible community - 1
set.seed(1651)

# creating 20 communities with the same Z1_range and z2_range
communities3 <- rerun(10, get_community(s, species_pars))

# add community ID  communities 1 

mynames <-  as.list(N_com)

communities3 <-map2(communities3, mynames, ~.x %>% mutate(community = .y))

(communities3_tibble <-  tibble(
  sp = map(communities3, "species"),
  E1 = map(communities3, "E1"),
  E2 = map(communities3, "E2"),
  rate = map(communities3, "rate"),
  community = map(communities3, "community")) %>%
    unnest(sp, E1, E2, rate, community))


these_E2 <- seq( min(communities3_tibble$E2), max(communities3_tibble$E2), length.out= 6)
temp <- communities3_tibble %>%
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

these_E1 <- seq( min(communities3_tibble$E1), max(communities3_tibble$E1), length.out= 6)
temp <- communities3_tibble %>%
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

performance_comm3 <- fig1 + fig2
performance_comm3

# get directional derivatives 
dir_3 <- map(communities3, ~ .x %>% 
               get_directionals(., new_data, refs))


mynames <-  as.list(N_com)

dir_3 <-map2(dir_3, mynames, ~.x %>% mutate(community = .y))

(dir_3.2 <- tibble(
  sp = map(dir_3, "sp"),
  time = map(dir_3, "time"),
  E1 = map(dir_3, "E1_ref"),
  E2 = map(dir_3, "E2_ref"),
  dir_deriv = map(dir_3, "dir_deriv"),
  community = map(dir_3, "community")) %>%
    unnest(sp, time, E1, E1, dir_deriv, community))

# from long to wide
rdiv_3 <- map(dir_3, ~ .x %>% 
                spread(key = sp, value = dir_deriv))


### Response diversity calculation
# dissimilarity
rdiv_3 <- lapply(rdiv_3, function(x) cbind(x,"rdiv"= apply(subset(x, select = c(s1, s2, s3, s4)), 1, resp_div, sign_sens = F)))
rdiv_3 <- lapply(rdiv_3, function(x) cbind(x,"Med"= apply(subset(x, select = c(rdiv)),2, median)))

# divergence
rdiv_3 <- lapply(rdiv_3, function(x) cbind(x,"sign"= apply(subset(x, select = c(s1, s2, s3, s4)), 1, resp_div, sign_sens = T)))
rdiv_3 <- lapply(rdiv_3, function(x) cbind(x,"Med_sing"= apply(subset(x, select = c(sign)),2, median)))


# from list to tibble
(rdiv_3.2 <- tibble(
  time = map(rdiv_3, "time"),
  E1 = map(rdiv_3, "E1_ref"),
  E2 = map(rdiv_3, "E2_ref"),
  rdiv = map(rdiv_3, "rdiv"),
  sign = map(rdiv_3, "sign"),
  community = map(rdiv_3, "community")) %>%
    unnest(time, E1, E2, rdiv, sign, community))

rdiv_3.2 <- rdiv_3.2 %>%  group_by(community) %>% 
  mutate(Med = mean(rdiv), 
         Med_sign = mean(sign))

# plotting 
rmax<-max(rdiv_3.2$rdiv); rmin<-1
dmax<-max(rdiv_3.2$sign); dmin<-0

Fig_1_comm3 <-dir_3.2 %>%
  ggplot(aes(x=time, y=dir_deriv, col = sp)) +
  #theme_classic(base_size = 14) + 
  labs(x = "time",y = "Directional derivative",tag = "a)") + 
  facet_wrap(vars(community), nrow = 2) +
  #geom_point(size=0.5) +
  geom_line() + theme_bw() + geom_hline(yintercept = 0, linetype= "dashed")


Fig2_comm3 <- rdiv_3.2 %>% 
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
                label = paste("RDiv^Mean =",paste0(round(rdiv_2.2$Med,digits = 4))))+
  lims(y = c(rmin,rmax))



Fig3_comm3 <- rdiv_3.2 %>% 
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
                label = paste("RDiv^Mean =",paste0(round(rdiv_2.2$Med_sign,digits = 4))))+
  lims(y = c(dmin,dmax))

Fig_1_comm3 / Fig2_comm3/ Fig3_comm3

### For summary plot

dd1 <- data.frame(rdiv_1.2[, c(4:6)]) 
dd1$cor <- "negative"
dd1$z1 <- "low"

dd2 <- data.frame(rdiv_2.2[, c(4:6)]) 
dd2$cor <- "negative"
dd2$z1 <- "medium"

dd3 <- data.frame(rdiv_3.2[, c(4:6)]) 
dd3$cor <- "negative"
dd3$z1 <- "high"




## High fluctuation in environmental change - positive correlation
#Increasing variance in z1 - positive correlation between env variables (E1 and E2) with *high* fluctuations

#Same steps as before (3 communities with increasing diversity ib z1 while z2 is fixed), but the correlation between E1 and E1 is negative.

library(MASS)
sample_size <- 15                                       
sample_medianvector <- c(mean(E1_series),  mean(E2_series))                                   
sample_covariance_matrix <- matrix(c(150, 150,
                                     150, 150),
                                   ncol = 2)

# create bivariate normal distribution
set.seed(5655)
refs <- mvrnorm(n = sample_size,
                mu = sample_meanvector, 
                Sigma = sample_covariance_matrix) %>%
  as_tibble() %>% 
  dplyr::rename(E1 = "V1", E2 = "V2") %>% 
  dplyr::mutate(time = seq.int(nrow(.)))



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



### Creating communities

# community 1 
s <- 4
# community 1 - low variation in z1. 
species_pars$z1_range <- 1.5
species_pars$z2_range <- 1.5

N_com <- seq(1, 10, 1)

# create reproducible community - 1
set.seed(1651)

# creating 20 communities with the same Z1_range and z2_range
communities1 <- rerun(10, get_community(s, species_pars))

# add community ID  communities 1 

mynames <-  as.list(N_com)

communities1 <-map2(communities1, mynames, ~.x %>% mutate(community = .y))
# get directional derivatives 
dir_1 <- map(communities1, ~ .x %>% 
               get_directionals(., new_data, refs))

mynames <-  as.list(N_com)

dir_1 <-map2(dir_1, mynames, ~.x %>% mutate(community = .y))


# community 2 
s <- 4
# community 1 - low variation in z1. 
species_pars$z1_range <- 25
species_pars$z2_range <- 25

N_com <- seq(1, 10, 1)

# create reproducible community - 2
set.seed(1651)

# creating 20 communities with the same Z1_range and z2_range
communities2 <- rerun(10, get_community(s, species_pars))

# add community ID  communities 2 

mynames <-  as.list(N_com)

communities2 <-map2(communities2, mynames, ~.x %>% mutate(community = .y))
# get directional derivatives 
dir_2 <- map(communities2, ~ .x %>% 
               get_directionals(., new_data, refs))


mynames <-  as.list(N_com)

dir_2 <-map2(dir_2, mynames, ~.x %>% mutate(community = .y))

# community 3
s <- 4
# community 1 - low variation in z1. 
species_pars$z1_range <- 50
species_pars$z2_range <- 50

N_com <- seq(1, 10, 1)

# create reproducible community - 1
set.seed(1651)

# creating 20 communities with the same Z1_range and z2_range
communities3 <- rerun(10, get_community(s, species_pars))

# add community ID  communities 1 

mynames <-  as.list(N_com)

communities3 <-map2(communities3, mynames, ~.x %>% mutate(community = .y))
# get directional derivatives 
dir_3 <- map(communities3, ~ .x %>% 
               get_directionals(., new_data, refs))


mynames <-  as.list(N_com)

dir_3 <-map2(dir_3, mynames, ~.x %>% mutate(community = .y))


### Response diversity calculation

#comm 1 
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



Fig3_comm1<- rdiv_1.2 %>% 
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


#comm 2
# from long to wide
rdiv_2 <- map(dir_2, ~ .x %>% 
                spread(key = sp, value = dir_deriv))


### Response diversity calculation
# dissimilarity
rdiv_2 <- lapply(rdiv_2, function(x) cbind(x,"rdiv"= apply(subset(x, select = c(s1, s2, s3, s4)), 1, resp_div, sign_sens = F)))
rdiv_2 <- lapply(rdiv_2, function(x) cbind(x,"Med"= apply(subset(x, select = c(rdiv)),2, median)))

# divergence
rdiv_2 <- lapply(rdiv_2, function(x) cbind(x,"sign"= apply(subset(x, select = c(s1, s2, s3, s4)), 1, resp_div, sign_sens = T)))
rdiv_2 <- lapply(rdiv_2, function(x) cbind(x,"Med_sing"= apply(subset(x, select = c(sign)),2, median)))


# from list to tibble
(rdiv_2.2 <- tibble(
  time = map(rdiv_2, "time"),
  E1 = map(rdiv_2, "E1_ref"),
  E2 = map(rdiv_2, "E2_ref"),
  rdiv = map(rdiv_2, "rdiv"),
  sign = map(rdiv_2, "sign"),
  community = map(rdiv_2, "community")) %>%
    unnest(time, E1, E2, rdiv, sign, community))

rdiv_2.2 <- rdiv_2.2 %>%  group_by(community) %>% 
  mutate(Med = mean(rdiv), 
         Med_sign = mean(sign))

# plotting 
rmax<-max(rdiv_2.2$rdiv); rmin<-1
dmax<-max(rdiv_2.2$sign); dmin<-0

Fig_1_comm2 <-dir_2.2 %>%
  ggplot(aes(x=time, y=dir_deriv, col = sp)) +
  #theme_classic(base_size = 14) + 
  labs(x = "time",y = "Directional derivative",tag = "a)") + 
  facet_wrap(vars(community), nrow = 2) +
  #geom_point(size=0.5) +
  geom_line() + theme_bw() + geom_hline(yintercept = 0, linetype= "dashed")


Fig2_comm2 <- rdiv_2.2 %>% 
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
                label = paste("RDiv^Mean =",paste0(round(rdiv_2.2$Med,digits = 4))))+
  lims(y = c(rmin,rmax))



Fig3_comm2<- rdiv_2.2 %>% 
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
                label = paste("RDiv^Mean =",paste0(round(rdiv_2.2$Med_sign,digits = 4))))+
  lims(y = c(dmin,dmax))

Fig_1_comm2 / Fig2_comm2/ Fig3_comm2




#comm 3
# from long to wide
rdiv_3 <- map(dir_3, ~ .x %>% 
                spread(key = sp, value = dir_deriv))


### Response diversity calculation
# dissimilarity
rdiv_3 <- lapply(rdiv_3, function(x) cbind(x,"rdiv"= apply(subset(x, select = c(s1, s2, s3, s4)), 1, resp_div, sign_sens = F)))
rdiv_3 <- lapply(rdiv_3, function(x) cbind(x,"Med"= apply(subset(x, select = c(rdiv)),2, median)))

# divergence
rdiv_3 <- lapply(rdiv_3, function(x) cbind(x,"sign"= apply(subset(x, select = c(s1, s2, s3, s4)), 1, resp_div, sign_sens = T)))
rdiv_3 <- lapply(rdiv_3, function(x) cbind(x,"Med_sing"= apply(subset(x, select = c(sign)),2, median)))


# from list to tibble
(rdiv_3.2 <- tibble(
  time = map(rdiv_3, "time"),
  E1 = map(rdiv_3, "E1_ref"),
  E2 = map(rdiv_3, "E2_ref"),
  rdiv = map(rdiv_3, "rdiv"),
  sign = map(rdiv_3, "sign"),
  community = map(rdiv_3, "community")) %>%
    unnest(time, E1, E2, rdiv, sign, community))

rdiv_3.2 <- rdiv_3.2 %>%  group_by(community) %>% 
  mutate(Med = mean(rdiv), 
         Med_sign = mean(sign))

# plotting 
rmax<-max(rdiv_3.2$rdiv); rmin<-1
dmax<-max(rdiv_3.2$sign); dmin<-0

Fig_1_comm3 <-dir_3.2 %>%
  ggplot(aes(x=time, y=dir_deriv, col = sp)) +
  #theme_classic(base_size = 14) + 
  labs(x = "time",y = "Directional derivative",tag = "a)") + 
  facet_wrap(vars(community), nrow = 2) +
  #geom_point(size=0.5) +
  geom_line() + theme_bw() + geom_hline(yintercept = 0, linetype= "dashed")


Fig2_comm3 <- rdiv_3.2 %>% 
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
                label = paste("RDiv^Mean =",paste0(round(rdiv_3.2$Med,digits = 4))))+
  lims(y = c(rmin,rmax))



Fig3_comm3<- rdiv_3.2 %>% 
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
                label = paste("RDiv^Mean =",paste0(round(rdiv_3.2$Med_sign,digits = 4))))+
  lims(y = c(dmin,dmax))

Fig_1_comm3/ Fig2_comm3/ Fig3_comm3


### For summary plot

dd4 <- data.frame(rdiv_1.2[, c(4:6)]) 
dd4$cor <- "positive"
dd4$z1 <- "low"

dd5 <- data.frame(rdiv_2.2[, c(4:6)]) 
dd5$cor <- "positive"
dd5$z1 <- "medium"

dd6 <- data.frame(rdiv_3.2[, c(4:6)]) 
dd6$cor <- "positive"
dd6$z1 <- "high"




## High fluctuation in environmental change - no correlation
# Increasing variance in z1 - no correlation between env variables (E1 and E2) with *high* fluctuations
# 
# Same steps as before (3 communities with increasing diversity in z1 while z2 is fixed), but there is no correlation between E1 and E1

library(MASS)
sample_size <- 15                                       
sample_meanvector <- c(mean(E1_series),  mean(E2_series))                                   
sample_covariance_matrix <- matrix(c(150, 0.0,
                                     0.0, 90),
                                   ncol = 2)

# create bivariate normal distribution
set.seed(5655)
refs <- mvrnorm(n = sample_size,
                mu = sample_meanvector, 
                Sigma = sample_covariance_matrix) %>%
  as_tibble() %>% 
  dplyr::rename(E1 = "V1", E2 = "V2") %>% 
  dplyr::mutate(time = seq.int(nrow(.)))


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






### Creating communities

# community 1 
s <- 4
# community 1 - low variation in z1. 
species_pars$z1_range <- 1.5
species_pars$z2_range <- 1.5

N_com <- seq(1, 10, 1)

# create reproducible community - 1
set.seed(1651)

# creating 20 communities with the same Z1_range and z2_range
communities1 <- rerun(10, get_community(s, species_pars))

# add community ID  communities 1 

mynames <-  as.list(N_com)

communities1 <-map2(communities1, mynames, ~.x %>% mutate(community = .y))
# get directional derivatives 
dir_1 <- map(communities1, ~ .x %>% 
               get_directionals(., new_data, refs))

mynames <-  as.list(N_com)

dir_1 <-map2(dir_1, mynames, ~.x %>% mutate(community = .y))


# community 2 
s <- 4
# community 1 - low variation in z1. 
species_pars$z1_range <- 25
species_pars$z2_range <- 25

N_com <- seq(1, 10, 1)

# create reproducible community - 2
set.seed(1651)

# creating 20 communities with the same Z1_range and z2_range
communities2 <- rerun(10, get_community(s, species_pars))

# add community ID  communities 2 

mynames <-  as.list(N_com)

communities2 <-map2(communities2, mynames, ~.x %>% mutate(community = .y))
# get directional derivatives 
dir_2 <- map(communities2, ~ .x %>% 
               get_directionals(., new_data, refs))


mynames <-  as.list(N_com)

dir_2 <-map2(dir_2, mynames, ~.x %>% mutate(community = .y))

# community 3
s <- 4
# community 1 - low variation in z1. 
species_pars$z1_range <- 50
species_pars$z2_range <- 50

N_com <- seq(1, 10, 1)

# create reproducible community - 1
set.seed(1651)

# creating 20 communities with the same Z1_range and z2_range
communities3 <- rerun(10, get_community(s, species_pars))

# add community ID  communities 1 

mynames <-  as.list(N_com)

communities3 <-map2(communities3, mynames, ~.x %>% mutate(community = .y))
# get directional derivatives 
dir_3 <- map(communities3, ~ .x %>% 
               get_directionals(., new_data, refs))


mynames <-  as.list(N_com)

dir_3 <-map2(dir_3, mynames, ~.x %>% mutate(community = .y))


### Response diversity calculation

#comm 1 
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



Fig3_comm1<- rdiv_1.2 %>% 
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


#comm 2
# from long to wide
rdiv_2 <- map(dir_2, ~ .x %>% 
                spread(key = sp, value = dir_deriv))


### Response diversity calculation
# dissimilarity
rdiv_2 <- lapply(rdiv_2, function(x) cbind(x,"rdiv"= apply(subset(x, select = c(s1, s2, s3, s4)), 1, resp_div, sign_sens = F)))
rdiv_2 <- lapply(rdiv_2, function(x) cbind(x,"Med"= apply(subset(x, select = c(rdiv)),2, median)))

# divergence
rdiv_2 <- lapply(rdiv_2, function(x) cbind(x,"sign"= apply(subset(x, select = c(s1, s2, s3, s4)), 1, resp_div, sign_sens = T)))
rdiv_2 <- lapply(rdiv_2, function(x) cbind(x,"Med_sing"= apply(subset(x, select = c(sign)),2, median)))


# from list to tibble
(rdiv_2.2 <- tibble(
  time = map(rdiv_2, "time"),
  E1 = map(rdiv_2, "E1_ref"),
  E2 = map(rdiv_2, "E2_ref"),
  rdiv = map(rdiv_2, "rdiv"),
  sign = map(rdiv_2, "sign"),
  community = map(rdiv_2, "community")) %>%
    unnest(time, E1, E2, rdiv, sign, community))

rdiv_2.2 <- rdiv_2.2 %>%  group_by(community) %>% 
  mutate(Med = mean(rdiv), 
         Med_sign = mean(sign))

# plotting 
rmax<-max(rdiv_2.2$rdiv); rmin<-1
dmax<-max(rdiv_2.2$sign); dmin<-0

Fig_1_comm2 <-dir_2.2 %>%
  ggplot(aes(x=time, y=dir_deriv, col = sp)) +
  #theme_classic(base_size = 14) + 
  labs(x = "time",y = "Directional derivative",tag = "a)") + 
  facet_wrap(vars(community), nrow = 2) +
  #geom_point(size=0.5) +
  geom_line() + theme_bw() + geom_hline(yintercept = 0, linetype= "dashed")


Fig2_comm2 <- rdiv_2.2 %>% 
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
                label = paste("RDiv^Mean =",paste0(round(rdiv_2.2$Med,digits = 4))))+
  lims(y = c(rmin,rmax))



Fig3_comm2<- rdiv_2.2 %>% 
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
                label = paste("RDiv^Mean =",paste0(round(rdiv_2.2$Med_sign,digits = 4))))+
  lims(y = c(dmin,dmax))

Fig_1_comm2 / Fig2_comm2/ Fig3_comm2




#comm 3
# from long to wide
rdiv_3 <- map(dir_3, ~ .x %>% 
                spread(key = sp, value = dir_deriv))


### Response diversity calculation
# dissimilarity
rdiv_3 <- lapply(rdiv_3, function(x) cbind(x,"rdiv"= apply(subset(x, select = c(s1, s2, s3, s4)), 1, resp_div, sign_sens = F)))
rdiv_3 <- lapply(rdiv_3, function(x) cbind(x,"Med"= apply(subset(x, select = c(rdiv)),2, median)))

# divergence
rdiv_3 <- lapply(rdiv_3, function(x) cbind(x,"sign"= apply(subset(x, select = c(s1, s2, s3, s4)), 1, resp_div, sign_sens = T)))
rdiv_3 <- lapply(rdiv_3, function(x) cbind(x,"Med_sing"= apply(subset(x, select = c(sign)),2, median)))


# from list to tibble
(rdiv_3.2 <- tibble(
  time = map(rdiv_3, "time"),
  E1 = map(rdiv_3, "E1_ref"),
  E2 = map(rdiv_3, "E2_ref"),
  rdiv = map(rdiv_3, "rdiv"),
  sign = map(rdiv_3, "sign"),
  community = map(rdiv_3, "community")) %>%
    unnest(time, E1, E2, rdiv, sign, community))

rdiv_3.2 <- rdiv_3.2 %>%  group_by(community) %>% 
  mutate(Med = mean(rdiv), 
         Med_sign = mean(sign))

# plotting 
rmax<-max(rdiv_3.2$rdiv); rmin<-1
dmax<-max(rdiv_3.2$sign); dmin<-0

Fig_1_comm3 <-dir_3.2 %>%
  ggplot(aes(x=time, y=dir_deriv, col = sp)) +
  #theme_classic(base_size = 14) + 
  labs(x = "time",y = "Directional derivative",tag = "a)") + 
  facet_wrap(vars(community), nrow = 2) +
  #geom_point(size=0.5) +
  geom_line() + theme_bw() + geom_hline(yintercept = 0, linetype= "dashed")


Fig2_comm3 <- rdiv_3.2 %>% 
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
                label = paste("RDiv^Mean =",paste0(round(rdiv_3.2$Med,digits = 4))))+
  lims(y = c(rmin,rmax))



Fig3_comm3<- rdiv_3.2 %>% 
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
                label = paste("RDiv^Mean =",paste0(round(rdiv_3.2$Med_sign,digits = 4))))+
  lims(y = c(dmin,dmax))

Fig_1_comm3/ Fig2_comm3/ Fig3_comm3


### For summary plot

dd7 <- data.frame(rdiv_1.2[, c(4:6)]) 
dd7$cor <- "no_cor"
dd7$z1 <- "low"

dd8 <- data.frame(rdiv_2.2[, c(4:6)]) 
dd8$cor <- "no_cor"
dd8$z1 <- "medium"

dd9 <- data.frame(rdiv_3.2[, c(4:6)]) 
dd9$cor <- "no_cor"
dd9$z1 <- "high"



dd_plot1 <- rbind(dd1, dd2, dd3,
                  dd4, dd5, dd6, 
                  dd7, dd8, dd9)

theme_set(theme_classic(base_size = 12))

g1.1 <-
  ggplot(dd_plot1, aes(x = cor, y = rdiv, color = z1)) +
  # scale_y_continuous(limits = c(0.99, 1.005), expand = c(0.02, 0.02)) +
  scale_color_uchicago() +
  labs(x = NULL, y = "Dissimilarity (derivatives)") +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 16),
    panel.grid = element_blank()
  )


rdiv2 <- g1.1 +
  geom_jitter(size = 2, alpha = 0.25, width = 0.2) +
  stat_summary(fun = mean, geom = "point", size = 5) + labs(tag = "(a)")



g1.2 <-
  ggplot(dd_plot1, aes(x = cor, y = sign, color = z1)) +
  scale_color_uchicago() +
  labs(x = NULL, y = "Divergence (sign sensitive)") +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 16),
    panel.grid = element_blank()
  )

sign2 <- g1.2 +
  geom_jitter(size = 2, alpha = 0.25, width = 0.2) +
  stat_summary(fun = mean, geom = "point", size = 5)+ labs(tag = "(b)") +
  ylim(0, 1)

combined <- rdiv1 + sign1 & theme(legend.position = "bottom")
combined + plot_layout(guides = "collect")











#######################################################
###### Variation in z1 and z2 - negative correlation ####
#######################################################
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
species_pars <- list(a1_mean = 1e-9,
                     b1_mean = 0.063,
                     z1_mean = 285,
                     w1_mean = 60,
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
species_pars$z2_range <- 50

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


these_E2 <- seq( min(communities1_tibble$E2), max(communities1_tibble$E2), length.out= 6)
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




#### Medium variation in z1 and z2 - negative correlation E1 E2 
s <- 4
# community 1 - low variation in z1. 
species_pars$z1_range <- 25
species_pars$z2_range <- 25

N_com <- seq(1, 10, 1)

# create reproducible community - 1
set.seed(1651)

# creating 20 communities with the same Z1_range and z2_range
communities2 <- rerun(10, get_community(s, species_pars))

# add community ID  communities 1 

mynames <-  as.list(N_com)

communities2 <-map2(communities2, mynames, ~.x %>% mutate(community = .y))

(communities2_tibble <-  tibble(
  sp = map(communities2, "species"),
  E1 = map(communities2, "E1"),
  E2 = map(communities2, "E2"),
  rate = map(communities2, "rate"),
  community = map(communities2, "community")) %>%
    unnest(sp, E1, E2, rate, community))


these_E2 <- seq( min(communities2_tibble$E2), max(communities2_tibble$E2), length.out= 6)
temp <- communities2_tibble %>%
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

these_E1 <- seq( min(communities2_tibble$E1), max(communities2_tibble$E1), length.out= 6)
temp <- communities2_tibble %>%
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

performance_comm2 <- fig1 + fig2
performance_comm2

# get directional derivatives 
dir_2 <- map(communities2, ~ .x %>% 
               get_directionals(., new_data, refs))


mynames <-  as.list(N_com)

dir_2 <-map2(dir_2, mynames, ~.x %>% mutate(community = .y))

(dir_2.2 <- tibble(
  sp = map(dir_2, "sp"),
  time = map(dir_2, "time"),
  E1 = map(dir_2, "E1_ref"),
  E2 = map(dir_2, "E2_ref"),
  dir_deriv = map(dir_2, "dir_deriv"),
  community = map(dir_2, "community")) %>%
    unnest(sp, time, E1, E1, dir_deriv, community))

# from long to wide
rdiv_2 <- map(dir_2, ~ .x %>% 
                spread(key = sp, value = dir_deriv))


### Response diversity calculation
# dissimilarity
rdiv_2 <- lapply(rdiv_2, function(x) cbind(x,"rdiv"= apply(subset(x, select = c(s1, s2, s3, s4)), 1, resp_div, sign_sens = F)))
rdiv_2 <- lapply(rdiv_2, function(x) cbind(x,"Med"= apply(subset(x, select = c(rdiv)),2, median)))

# divergence
rdiv_2 <- lapply(rdiv_2, function(x) cbind(x,"sign"= apply(subset(x, select = c(s1, s2, s3, s4)), 1, resp_div, sign_sens = T)))
rdiv_2 <- lapply(rdiv_2, function(x) cbind(x,"Med_sing"= apply(subset(x, select = c(sign)),2, median)))


# from list to tibble
(rdiv_2.2 <- tibble(
  time = map(rdiv_2, "time"),
  E1 = map(rdiv_2, "E1_ref"),
  E2 = map(rdiv_2, "E2_ref"),
  rdiv = map(rdiv_2, "rdiv"),
  sign = map(rdiv_2, "sign"),
  community = map(rdiv_2, "community")) %>%
    unnest(time, E1, E2, rdiv, sign, community))

rdiv_2.2 <- rdiv_2.2 %>%  group_by(community) %>% 
  mutate(Med = mean(rdiv), 
         Med_sign = mean(sign))

# plotting 
rmax<-max(rdiv_2.2$rdiv); rmin<-1
dmax<-max(rdiv_2.2$sign); dmin<-0

Fig_1_comm2 <-dir_2.2 %>%
  ggplot(aes(x=time, y=dir_deriv, col = sp)) +
  #theme_classic(base_size = 14) + 
  labs(x = "time",y = "Directional derivative",tag = "a)") + 
  facet_wrap(vars(community), nrow = 2) +
  #geom_point(size=0.5) +
  geom_line() + theme_bw() + geom_hline(yintercept = 0, linetype= "dashed")


Fig2_comm2 <- rdiv_2.2 %>% 
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
                label = paste("RDiv^Mean =",paste0(round(rdiv_2.2$Med,digits = 4))))+
  lims(y = c(rmin,rmax))



Fig3_comm2 <- rdiv_2.2 %>% 
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
                label = paste("RDiv^Mean =",paste0(round(rdiv_2.2$Med_sign,digits = 4))))+
  lims(y = c(dmin,dmax))

Fig_1_comm2 / Fig2_comm2/ Fig3_comm2





#### High variation in z1 - negative correlation E1 E2 
s <- 4
# community 1 - low variation in z1. 
species_pars$z1_range <- 50
species_pars$z2_range <- 1.5

N_com <- seq(1, 10, 1)

# create reproducible community - 1
set.seed(1651)

# creating 20 communities with the same Z1_range and z2_range
communities3 <- rerun(10, get_community(s, species_pars))

# add community ID  communities 1 

mynames <-  as.list(N_com)

communities3 <-map2(communities3, mynames, ~.x %>% mutate(community = .y))

(communities3_tibble <-  tibble(
  sp = map(communities3, "species"),
  E1 = map(communities3, "E1"),
  E2 = map(communities3, "E2"),
  rate = map(communities3, "rate"),
  community = map(communities3, "community")) %>%
    unnest(sp, E1, E2, rate, community))


these_E2 <- seq( min(communities3_tibble$E2), max(communities3_tibble$E2), length.out= 6)
temp <- communities3_tibble %>%
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

these_E1 <- seq( min(communities3_tibble$E1), max(communities3_tibble$E1), length.out= 6)
temp <- communities3_tibble %>%
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

performance_comm3 <- fig1 + fig2
performance_comm3

# get directional derivatives 
dir_3 <- map(communities3, ~ .x %>% 
               get_directionals(., new_data, refs))


mynames <-  as.list(N_com)

dir_3 <-map2(dir_3, mynames, ~.x %>% mutate(community = .y))

(dir_3.2 <- tibble(
  sp = map(dir_3, "sp"),
  time = map(dir_3, "time"),
  E1 = map(dir_3, "E1_ref"),
  E2 = map(dir_3, "E2_ref"),
  dir_deriv = map(dir_3, "dir_deriv"),
  community = map(dir_3, "community")) %>%
    unnest(sp, time, E1, E1, dir_deriv, community))

# from long to wide
rdiv_3 <- map(dir_3, ~ .x %>% 
                spread(key = sp, value = dir_deriv))


### Response diversity calculation
# dissimilarity
rdiv_3 <- lapply(rdiv_3, function(x) cbind(x,"rdiv"= apply(subset(x, select = c(s1, s2, s3, s4)), 1, resp_div, sign_sens = F)))
rdiv_3 <- lapply(rdiv_3, function(x) cbind(x,"Med"= apply(subset(x, select = c(rdiv)),2, median)))

# divergence
rdiv_3 <- lapply(rdiv_3, function(x) cbind(x,"sign"= apply(subset(x, select = c(s1, s2, s3, s4)), 1, resp_div, sign_sens = T)))
rdiv_3 <- lapply(rdiv_3, function(x) cbind(x,"Med_sing"= apply(subset(x, select = c(sign)),2, median)))


# from list to tibble
(rdiv_3.2 <- tibble(
  time = map(rdiv_3, "time"),
  E1 = map(rdiv_3, "E1_ref"),
  E2 = map(rdiv_3, "E2_ref"),
  rdiv = map(rdiv_3, "rdiv"),
  sign = map(rdiv_3, "sign"),
  community = map(rdiv_3, "community")) %>%
    unnest(time, E1, E2, rdiv, sign, community))

rdiv_3.2 <- rdiv_3.2 %>%  group_by(community) %>% 
  mutate(Med = mean(rdiv), 
         Med_sign = mean(sign))

# plotting 
rmax<-max(rdiv_3.2$rdiv); rmin<-1
dmax<-max(rdiv_3.2$sign); dmin<-0

Fig_1_comm3 <-dir_3.2 %>%
  ggplot(aes(x=time, y=dir_deriv, col = sp)) +
  #theme_classic(base_size = 14) + 
  labs(x = "time",y = "Directional derivative",tag = "a)") + 
  facet_wrap(vars(community), nrow = 2) +
  #geom_point(size=0.5) +
  geom_line() + theme_bw() + geom_hline(yintercept = 0, linetype= "dashed")


Fig2_comm3 <- rdiv_3.2 %>% 
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
                label = paste("RDiv^Mean =",paste0(round(rdiv_2.2$Med,digits = 4))))+
  lims(y = c(rmin,rmax))



Fig3_comm3 <- rdiv_3.2 %>% 
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
                label = paste("RDiv^Mean =",paste0(round(rdiv_2.2$Med_sign,digits = 4))))+
  lims(y = c(dmin,dmax))

Fig_1_comm3 / Fig2_comm3/ Fig3_comm3

### For summary plot

dd1 <- data.frame(rdiv_1.2[, c(4:6)]) 
dd1$cor <- "negative"
dd1$z1 <- "low"

dd2 <- data.frame(rdiv_2.2[, c(4:6)]) 
dd2$cor <- "negative"
dd2$z1 <- "medium"

dd3 <- data.frame(rdiv_3.2[, c(4:6)]) 
dd3$cor <- "negative"
dd3$z1 <- "high"




## High fluctuation in environmental change - positive correlation
#Increasing variance in z1 - positive correlation between env variables (E1 and E2) with *high* fluctuations

#Same steps as before (3 communities with increasing diversity ib z1 while z2 is fixed), but the correlation between E1 and E1 is negative.

library(MASS)
sample_size <- 15                                       
sample_medianvector <- c(mean(E1_series),  mean(E2_series))                                   
sample_covariance_matrix <- matrix(c(150, 150,
                                     150, 150),
                                   ncol = 2)

# create bivariate normal distribution
set.seed(5655)
refs <- mvrnorm(n = sample_size,
                mu = sample_meanvector, 
                Sigma = sample_covariance_matrix) %>%
  as_tibble() %>% 
  dplyr::rename(E1 = "V1", E2 = "V2") %>% 
  dplyr::mutate(time = seq.int(nrow(.)))



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



### Creating communities

# community 1 
s <- 4
# community 1 - low variation in z1. 
species_pars$z1_range <- 1.5
species_pars$z2_range <- 50

N_com <- seq(1, 10, 1)

# create reproducible community - 1
set.seed(1651)

# creating 20 communities with the same Z1_range and z2_range
communities1 <- rerun(10, get_community(s, species_pars))

# add community ID  communities 1 

mynames <-  as.list(N_com)

communities1 <-map2(communities1, mynames, ~.x %>% mutate(community = .y))
# get directional derivatives 
dir_1 <- map(communities1, ~ .x %>% 
               get_directionals(., new_data, refs))

mynames <-  as.list(N_com)

dir_1 <-map2(dir_1, mynames, ~.x %>% mutate(community = .y))


# community 2 
s <- 4
# community 1 - low variation in z1. 
species_pars$z1_range <- 25
species_pars$z2_range <- 25

N_com <- seq(1, 10, 1)

# create reproducible community - 2
set.seed(1651)

# creating 20 communities with the same Z1_range and z2_range
communities2 <- rerun(10, get_community(s, species_pars))

# add community ID  communities 2 

mynames <-  as.list(N_com)

communities2 <-map2(communities2, mynames, ~.x %>% mutate(community = .y))
# get directional derivatives 
dir_2 <- map(communities2, ~ .x %>% 
               get_directionals(., new_data, refs))


mynames <-  as.list(N_com)

dir_2 <-map2(dir_2, mynames, ~.x %>% mutate(community = .y))

# community 3
s <- 4
# community 1 - low variation in z1. 
species_pars$z1_range <- 50
species_pars$z2_range <- 1.5

N_com <- seq(1, 10, 1)

# create reproducible community - 1
set.seed(1651)

# creating 20 communities with the same Z1_range and z2_range
communities3 <- rerun(10, get_community(s, species_pars))

# add community ID  communities 1 

mynames <-  as.list(N_com)

communities3 <-map2(communities3, mynames, ~.x %>% mutate(community = .y))
# get directional derivatives 
dir_3 <- map(communities3, ~ .x %>% 
               get_directionals(., new_data, refs))


mynames <-  as.list(N_com)

dir_3 <-map2(dir_3, mynames, ~.x %>% mutate(community = .y))


### Response diversity calculation

#comm 1 
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



Fig3_comm1<- rdiv_1.2 %>% 
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


#comm 2
# from long to wide
rdiv_2 <- map(dir_2, ~ .x %>% 
                spread(key = sp, value = dir_deriv))


### Response diversity calculation
# dissimilarity
rdiv_2 <- lapply(rdiv_2, function(x) cbind(x,"rdiv"= apply(subset(x, select = c(s1, s2, s3, s4)), 1, resp_div, sign_sens = F)))
rdiv_2 <- lapply(rdiv_2, function(x) cbind(x,"Med"= apply(subset(x, select = c(rdiv)),2, median)))

# divergence
rdiv_2 <- lapply(rdiv_2, function(x) cbind(x,"sign"= apply(subset(x, select = c(s1, s2, s3, s4)), 1, resp_div, sign_sens = T)))
rdiv_2 <- lapply(rdiv_2, function(x) cbind(x,"Med_sing"= apply(subset(x, select = c(sign)),2, median)))


# from list to tibble
(rdiv_2.2 <- tibble(
  time = map(rdiv_2, "time"),
  E1 = map(rdiv_2, "E1_ref"),
  E2 = map(rdiv_2, "E2_ref"),
  rdiv = map(rdiv_2, "rdiv"),
  sign = map(rdiv_2, "sign"),
  community = map(rdiv_2, "community")) %>%
    unnest(time, E1, E2, rdiv, sign, community))

rdiv_2.2 <- rdiv_2.2 %>%  group_by(community) %>% 
  mutate(Med = mean(rdiv), 
         Med_sign = mean(sign))

# plotting 
rmax<-max(rdiv_2.2$rdiv); rmin<-1
dmax<-max(rdiv_2.2$sign); dmin<-0

Fig_1_comm2 <-dir_2.2 %>%
  ggplot(aes(x=time, y=dir_deriv, col = sp)) +
  #theme_classic(base_size = 14) + 
  labs(x = "time",y = "Directional derivative",tag = "a)") + 
  facet_wrap(vars(community), nrow = 2) +
  #geom_point(size=0.5) +
  geom_line() + theme_bw() + geom_hline(yintercept = 0, linetype= "dashed")


Fig2_comm2 <- rdiv_2.2 %>% 
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
                label = paste("RDiv^Mean =",paste0(round(rdiv_2.2$Med,digits = 4))))+
  lims(y = c(rmin,rmax))



Fig3_comm2<- rdiv_2.2 %>% 
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
                label = paste("RDiv^Mean =",paste0(round(rdiv_2.2$Med_sign,digits = 4))))+
  lims(y = c(dmin,dmax))

Fig_1_comm2 / Fig2_comm2/ Fig3_comm2




#comm 3
# from long to wide
rdiv_3 <- map(dir_3, ~ .x %>% 
                spread(key = sp, value = dir_deriv))


### Response diversity calculation
# dissimilarity
rdiv_3 <- lapply(rdiv_3, function(x) cbind(x,"rdiv"= apply(subset(x, select = c(s1, s2, s3, s4)), 1, resp_div, sign_sens = F)))
rdiv_3 <- lapply(rdiv_3, function(x) cbind(x,"Med"= apply(subset(x, select = c(rdiv)),2, median)))

# divergence
rdiv_3 <- lapply(rdiv_3, function(x) cbind(x,"sign"= apply(subset(x, select = c(s1, s2, s3, s4)), 1, resp_div, sign_sens = T)))
rdiv_3 <- lapply(rdiv_3, function(x) cbind(x,"Med_sing"= apply(subset(x, select = c(sign)),2, median)))


# from list to tibble
(rdiv_3.2 <- tibble(
  time = map(rdiv_3, "time"),
  E1 = map(rdiv_3, "E1_ref"),
  E2 = map(rdiv_3, "E2_ref"),
  rdiv = map(rdiv_3, "rdiv"),
  sign = map(rdiv_3, "sign"),
  community = map(rdiv_3, "community")) %>%
    unnest(time, E1, E2, rdiv, sign, community))

rdiv_3.2 <- rdiv_3.2 %>%  group_by(community) %>% 
  mutate(Med = mean(rdiv), 
         Med_sign = mean(sign))

# plotting 
rmax<-max(rdiv_3.2$rdiv); rmin<-1
dmax<-max(rdiv_3.2$sign); dmin<-0

Fig_1_comm3 <-dir_3.2 %>%
  ggplot(aes(x=time, y=dir_deriv, col = sp)) +
  #theme_classic(base_size = 14) + 
  labs(x = "time",y = "Directional derivative",tag = "a)") + 
  facet_wrap(vars(community), nrow = 2) +
  #geom_point(size=0.5) +
  geom_line() + theme_bw() + geom_hline(yintercept = 0, linetype= "dashed")


Fig2_comm3 <- rdiv_3.2 %>% 
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
                label = paste("RDiv^Mean =",paste0(round(rdiv_3.2$Med,digits = 4))))+
  lims(y = c(rmin,rmax))



Fig3_comm3<- rdiv_3.2 %>% 
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
                label = paste("RDiv^Mean =",paste0(round(rdiv_3.2$Med_sign,digits = 4))))+
  lims(y = c(dmin,dmax))

Fig_1_comm3/ Fig2_comm3/ Fig3_comm3


### For summary plot

dd4 <- data.frame(rdiv_1.2[, c(4:6)]) 
dd4$cor <- "positive"
dd4$z1 <- "low"

dd5 <- data.frame(rdiv_2.2[, c(4:6)]) 
dd5$cor <- "positive"
dd5$z1 <- "medium"

dd6 <- data.frame(rdiv_3.2[, c(4:6)]) 
dd6$cor <- "positive"
dd6$z1 <- "high"




## High fluctuation in environmental change - no correlation
# Increasing variance in z1 - no correlation between env variables (E1 and E2) with *high* fluctuations
# 
# Same steps as before (3 communities with increasing diversity in z1 while z2 is fixed), but there is no correlation between E1 and E1

library(MASS)
sample_size <- 15                                       
sample_meanvector <- c(mean(E1_series),  mean(E2_series))                                   
sample_covariance_matrix <- matrix(c(150, 0.0,
                                     0.0, 90),
                                   ncol = 2)

# create bivariate normal distribution
set.seed(5655)
refs <- mvrnorm(n = sample_size,
                mu = sample_meanvector, 
                Sigma = sample_covariance_matrix) %>%
  as_tibble() %>% 
  dplyr::rename(E1 = "V1", E2 = "V2") %>% 
  dplyr::mutate(time = seq.int(nrow(.)))


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






### Creating communities

# community 1 
s <- 4
# community 1 - low variation in z1. 
species_pars$z1_range <- 1.5
species_pars$z2_range <- 50

N_com <- seq(1, 10, 1)

# create reproducible community - 1
set.seed(1651)

# creating 20 communities with the same Z1_range and z2_range
communities1 <- rerun(10, get_community(s, species_pars))

# add community ID  communities 1 

mynames <-  as.list(N_com)

communities1 <-map2(communities1, mynames, ~.x %>% mutate(community = .y))
# get directional derivatives 
dir_1 <- map(communities1, ~ .x %>% 
               get_directionals(., new_data, refs))

mynames <-  as.list(N_com)

dir_1 <-map2(dir_1, mynames, ~.x %>% mutate(community = .y))


# community 2 
s <- 4
# community 1 - low variation in z1. 
species_pars$z1_range <- 25
species_pars$z2_range <- 25

N_com <- seq(1, 10, 1)

# create reproducible community - 2
set.seed(1651)

# creating 20 communities with the same Z1_range and z2_range
communities2 <- rerun(10, get_community(s, species_pars))

# add community ID  communities 2 

mynames <-  as.list(N_com)

communities2 <-map2(communities2, mynames, ~.x %>% mutate(community = .y))
# get directional derivatives 
dir_2 <- map(communities2, ~ .x %>% 
               get_directionals(., new_data, refs))


mynames <-  as.list(N_com)

dir_2 <-map2(dir_2, mynames, ~.x %>% mutate(community = .y))

# community 3
s <- 4
# community 1 - low variation in z1. 
species_pars$z1_range <- 50
species_pars$z2_range <- 1.5

N_com <- seq(1, 10, 1)

# create reproducible community - 1
set.seed(1651)

# creating 20 communities with the same Z1_range and z2_range
communities3 <- rerun(10, get_community(s, species_pars))

# add community ID  communities 1 

mynames <-  as.list(N_com)

communities3 <-map2(communities3, mynames, ~.x %>% mutate(community = .y))
# get directional derivatives 
dir_3 <- map(communities3, ~ .x %>% 
               get_directionals(., new_data, refs))


mynames <-  as.list(N_com)

dir_3 <-map2(dir_3, mynames, ~.x %>% mutate(community = .y))


### Response diversity calculation

#comm 1 
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



Fig3_comm1<- rdiv_1.2 %>% 
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


#comm 2
# from long to wide
rdiv_2 <- map(dir_2, ~ .x %>% 
                spread(key = sp, value = dir_deriv))


### Response diversity calculation
# dissimilarity
rdiv_2 <- lapply(rdiv_2, function(x) cbind(x,"rdiv"= apply(subset(x, select = c(s1, s2, s3, s4)), 1, resp_div, sign_sens = F)))
rdiv_2 <- lapply(rdiv_2, function(x) cbind(x,"Med"= apply(subset(x, select = c(rdiv)),2, median)))

# divergence
rdiv_2 <- lapply(rdiv_2, function(x) cbind(x,"sign"= apply(subset(x, select = c(s1, s2, s3, s4)), 1, resp_div, sign_sens = T)))
rdiv_2 <- lapply(rdiv_2, function(x) cbind(x,"Med_sing"= apply(subset(x, select = c(sign)),2, median)))


# from list to tibble
(rdiv_2.2 <- tibble(
  time = map(rdiv_2, "time"),
  E1 = map(rdiv_2, "E1_ref"),
  E2 = map(rdiv_2, "E2_ref"),
  rdiv = map(rdiv_2, "rdiv"),
  sign = map(rdiv_2, "sign"),
  community = map(rdiv_2, "community")) %>%
    unnest(time, E1, E2, rdiv, sign, community))

rdiv_2.2 <- rdiv_2.2 %>%  group_by(community) %>% 
  mutate(Med = mean(rdiv), 
         Med_sign = mean(sign))

# plotting 
rmax<-max(rdiv_2.2$rdiv); rmin<-1
dmax<-max(rdiv_2.2$sign); dmin<-0

Fig_1_comm2 <-dir_2.2 %>%
  ggplot(aes(x=time, y=dir_deriv, col = sp)) +
  #theme_classic(base_size = 14) + 
  labs(x = "time",y = "Directional derivative",tag = "a)") + 
  facet_wrap(vars(community), nrow = 2) +
  #geom_point(size=0.5) +
  geom_line() + theme_bw() + geom_hline(yintercept = 0, linetype= "dashed")


Fig2_comm2 <- rdiv_2.2 %>% 
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
                label = paste("RDiv^Mean =",paste0(round(rdiv_2.2$Med,digits = 4))))+
  lims(y = c(rmin,rmax))



Fig3_comm2<- rdiv_2.2 %>% 
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
                label = paste("RDiv^Mean =",paste0(round(rdiv_2.2$Med_sign,digits = 4))))+
  lims(y = c(dmin,dmax))

Fig_1_comm2 / Fig2_comm2/ Fig3_comm2




#comm 3
# from long to wide
rdiv_3 <- map(dir_3, ~ .x %>% 
                spread(key = sp, value = dir_deriv))


### Response diversity calculation
# dissimilarity
rdiv_3 <- lapply(rdiv_3, function(x) cbind(x,"rdiv"= apply(subset(x, select = c(s1, s2, s3, s4)), 1, resp_div, sign_sens = F)))
rdiv_3 <- lapply(rdiv_3, function(x) cbind(x,"Med"= apply(subset(x, select = c(rdiv)),2, median)))

# divergence
rdiv_3 <- lapply(rdiv_3, function(x) cbind(x,"sign"= apply(subset(x, select = c(s1, s2, s3, s4)), 1, resp_div, sign_sens = T)))
rdiv_3 <- lapply(rdiv_3, function(x) cbind(x,"Med_sing"= apply(subset(x, select = c(sign)),2, median)))


# from list to tibble
(rdiv_3.2 <- tibble(
  time = map(rdiv_3, "time"),
  E1 = map(rdiv_3, "E1_ref"),
  E2 = map(rdiv_3, "E2_ref"),
  rdiv = map(rdiv_3, "rdiv"),
  sign = map(rdiv_3, "sign"),
  community = map(rdiv_3, "community")) %>%
    unnest(time, E1, E2, rdiv, sign, community))

rdiv_3.2 <- rdiv_3.2 %>%  group_by(community) %>% 
  mutate(Med = mean(rdiv), 
         Med_sign = mean(sign))

# plotting 
rmax<-max(rdiv_3.2$rdiv); rmin<-1
dmax<-max(rdiv_3.2$sign); dmin<-0

Fig_1_comm3 <-dir_3.2 %>%
  ggplot(aes(x=time, y=dir_deriv, col = sp)) +
  #theme_classic(base_size = 14) + 
  labs(x = "time",y = "Directional derivative",tag = "a)") + 
  facet_wrap(vars(community), nrow = 2) +
  #geom_point(size=0.5) +
  geom_line() + theme_bw() + geom_hline(yintercept = 0, linetype= "dashed")


Fig2_comm3 <- rdiv_3.2 %>% 
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
                label = paste("RDiv^Mean =",paste0(round(rdiv_3.2$Med,digits = 4))))+
  lims(y = c(rmin,rmax))



Fig3_comm3<- rdiv_3.2 %>% 
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
                label = paste("RDiv^Mean =",paste0(round(rdiv_3.2$Med_sign,digits = 4))))+
  lims(y = c(dmin,dmax))

Fig_1_comm3/ Fig2_comm3/ Fig3_comm3


### For summary plot

dd7 <- data.frame(rdiv_1.2[, c(4:6)]) 
dd7$cor <- "no_cor"
dd7$z1 <- "low"

dd8 <- data.frame(rdiv_2.2[, c(4:6)]) 
dd8$cor <- "no_cor"
dd8$z1 <- "medium"

dd9 <- data.frame(rdiv_3.2[, c(4:6)]) 
dd9$cor <- "no_cor"
dd9$z1 <- "high"



dd_plot1 <- rbind(dd1, dd2, dd3,
                  dd4, dd5, dd6, 
                  dd7, dd8, dd9)

theme_set(theme_classic(base_size = 12))

g1.1 <-
  ggplot(dd_plot1, aes(x = cor, y = rdiv, color = z1)) +
  # scale_y_continuous(limits = c(0.99, 1.005), expand = c(0.02, 0.02)) +
  scale_color_uchicago() +
  labs(x = NULL, y = "Dissimilarity (derivatives)") +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 16),
    panel.grid = element_blank()
  )


rdiv3 <- g1.1 +
  geom_jitter(size = 2, alpha = 0.25, width = 0.2) +
  stat_summary(fun = mean, geom = "point", size = 5) + labs(tag = "(a)")



g1.2 <-
  ggplot(dd_plot1, aes(x = cor, y = sign, color = z1)) +
  scale_color_uchicago() +
  labs(x = NULL, y = "Divergence (sign sensitive)") +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 16),
    panel.grid = element_blank()
  )

sign3 <- g1.2 +
  geom_jitter(size = 2, alpha = 0.25, width = 0.2) +
  stat_summary(fun = mean, geom = "point", size = 5)+ labs(tag = "(b)") +
  ylim(0, 1)

combined <- rdiv3 + sign3 & theme(legend.position = "bottom")
combined + plot_layout(guides = "collect")



rdiv1 <- rdiv1 + ggtitle("Only variation in z1") + labs(tag = "(a)")
rdiv2 <- rdiv2 + ggtitle("Variation in z1 and z2\npositive correlation") + labs(tag = "(b)")
rdiv3 <- rdiv3 + ggtitle("Variation in z1 and z2\nnegative correlation") + labs(tag = "(c)")

sign1 <- sign1 + ggtitle("Only variation in z1") + labs(tag = "(d)")
sign2 <- sign2 + ggtitle("Variation in z1 and z2\npositive correlation") + labs(tag = "(e)")
sign3 <- sign3 + ggtitle("Variation in z1 and z2\nnegative correlation") + labs(tag = "(e)")

combined <- rdiv1 + rdiv2 + rdiv3 + sign1 + sign2 +sign3 & theme(legend.position = "bottom")
combined + plot_layout(guides = "collect")