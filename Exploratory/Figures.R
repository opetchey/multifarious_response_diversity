radius <- 1
num_arrows <- 3
(pd_spp_no <- tibble(angle = rep(seq(0, 2*pi, length = num_arrows), nrow(pd_spp_no))) %>%
  mutate(E1_ref = rep(pd_spp_no$E1_ref, each = num_arrows),
         E2_ref = rep(pd_spp_no$E2_ref, each = num_arrows)) %>%
  full_join(pd_spp_no) %>%
  mutate(E1 = cos(angle) * radius,
         E2 = sin(angle) * radius,
         dir_deriv = E1 * pd_E1 + E2 * pd_E2,
         unit_vec_mag = sqrt(E1^2 + E2^2)))


# filter useless stuff
(pd_spp_no <- pd_spp_no %>% 
  dplyr::select(angle, sp, E1_ref, E2_ref, dir_deriv))


pd_spp_no %>% dplyr::filter(sp == "s1")
# Convert df into a list based on angle (so direction of env change)
(lst_dd <- pd_spp_no %>% group_by(angle) %>% 
  group_map(~.x))

# add column with point id to identify points on the surface
point_len <- 4
(lst_dd <- map(lst_dd, ~ .x %>% 
                mutate(point = rep(seq(1, 1 + nrow(.x) %/% point_len), each = point_len, length.out = nrow(.x)))))


(lst_dd.2 <- lst_dd %>% map(~ .x %>% 
                             dplyr::group_by(sp) %>% 
                             summarise_at(vars(-point), list(mean = mean)) %>% 
                             select(-c(E1_ref_mean, E2_ref_mean))))


try_1 <- tibble(
  sp = map(lst_dd.2, "sp"),
  dir_deriv = map(lst_dd.2, "dir_deriv_mean")) %>% 
  mutate(angle = row_number()) %>% 
  unnest(sp, dir_deriv)

try_2 <- try_1 %>% dplyr::filter(sp == "s4")

try_2$E1_ref <- 300
try_2$E2_ref <- 20


dd1_try2<- try_2 %>% 
  mutate(E1 = cos(angle) * radius,
         E2 = sin(angle) * radius)

sp2 +  geom_segment(data = dd1_try2,
                    aes(x = E1_ref, y = E2_ref,
                        xend = E1_ref + E1, yend = E2_ref + E2,
                        col = dir_deriv), alpha = 0.7,
                    size=1) +  scale_colour_gradient2(low = muted("blue"),
                                                      mid = "white",
                                                      high = muted("red"),
                                                      midpoint = 0)+
  theme(legend.position = "none")

# Creating the dataset with the 2 columns as described above
predicted <- nested_gams_no_inter %>% unnest(predicted)
rates <- cbind(new_data, predicted[,c(1,4)])
rates <- rates %>%
  relocate(species, E1, E2, predicted)




avg <- rates %>%  dplyr::filter(species == "s1")



data <- data.frame(x=avg$E1,y=avg$E2,z=avg$predicted)

density_avg <- as.matrix(acast(data, x~y, value.var="z"))

myxticks = c(0.55, 1.0, 1.5, 2, 2.5, 3.0, 3.5, 4.0, 4.5, 5)
axx <- list(
  title = 'Temperature',
  nticks = myxticks,
  ticktext = myxticks
)

myyticks = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)
axy <- list(
  title = 'Salinity',
  nticks = myyticks,
  ticktext = myyticks
)

axz <- list(
  title = 'Growth rate',
  nticks = 10
)



x_Avg <- as.numeric(rownames(density_avg))
y_Avg <- as.numeric(colnames(density_avg))


figAvg <- plot_ly(z = ~density_avg) %>% 
  add_surface(contours = list(z = list(show=TRUE, usecolormap=TRUE, project=list(z=TRUE))))%>%
  layout(title = 'Growth rate sp 1', plot_bgcolor = "#e5ecf6", scene = list(xaxis=axx,yaxis=axy,zaxis=axz))

figAvg 




x <- c(1, 2, 3)
y <- c(4, 5, 6)
z <- c(7, 8, 9)

fig <- plot_ly() %>%
  add_trace(x = x, y = y, z = z, type = "scatter3d", mode = "lines") %>%
  add_trace(x = c(1, 1), y = c(4, 5), z = c(7, 8), type = "scatter3d", mode = "lines") %>%
  add_trace(x = c(1, 1), y = c(4, 4), z = c(7, 9), type = "scatter3d", mode = "lines") %>%
  add_trace(x = c(1, 2), y = c(4, 4), z = c(7, 7), type = "scatter3d", mode = "lines") %>%
  layout(scene = list(aspectmode = "data"))

fig


(final_lst <- lst_dd.2 %>%  map( ~ .x %>% 
                                  spread(key = sp, value = dir_deriv_mean)))


### Response diversity calculation
# dissimilarity
rdiv_1 <- lapply(final_lst, function(x) cbind(x,"rdiv"= apply(x, 1, resp_div, sign_sens = F)))


# divergence
rdiv_1 <- lapply(rdiv_1, function(x) cbind(x,"sign"= apply(x, 1, resp_div, sign_sens = T)))


# from list to tibble
rdiv_1 <- tibble(
  rdiv = map(rdiv_1, "rdiv"),
  sing = map(rdiv_1, "sign")) %>% 
  mutate(points = row_number()) %>% 
  unnest(rdiv, sing)

ggarrange(sp1, sp2, sp3, sp4, ncol=4, nrow=1, common.legend = TRUE, legend="right")







# Figure directional derivatives in one directions

refs <- crossing(E1 = 273.15 + seq(0, 50, 5),
                 E2 = seq(0, 50, 5))



m_list_no <- (nested_gams_no_inter$gams)

my_spp_names_no <- (nested_gams_no_inter$species)


(pd_list_no <- modify_depth(m_list_no, 1, ~ get_partials(., refs)))

(refs <- tibble(
  E1_ref = map(pd_list_no, "E1"),
  E2_ref = map(pd_list_no, "E2"),
  pd_E1 = map(pd_list_no, "pd_E1"),
  pd_E2 = map(pd_list_no, "pd_E2")) %>% 
    dplyr::mutate(sp = my_spp_names_no) %>% 
    relocate(sp, E1_ref, E2_ref, pd_E1, pd_E2) %>% 
    unnest(E1_ref, E2_ref, pd_E1, pd_E2))


radius <- 5
num_arrows <-1


refs_sp1 <- refs %>% dplyr::filter(sp == "s1")
dd1_sp1<- tibble(angle = rep(seq(0.8, 2*pi, length = num_arrows), nrow(refs_sp1))) %>% 
  mutate(E1_ref = rep(refs_sp1$E1_ref, each = num_arrows),
         E2_ref = rep(refs_sp1$E2_ref, each = num_arrows)) %>%  
  full_join(refs_sp1) %>%  
  mutate(E1 = cos(angle) * radius,
         E2 = sin(angle) * radius,
         dir_deriv = E1 * pd_E1 + E2 * pd_E2, unit_vec_mag = sqrt(E1^2 + E2^2))

sp1 <- sp1 +  geom_segment(data = dd1_sp1,
                           aes(x = E1_ref, y = E2_ref,
                               xend = E1_ref + E1, yend = E2_ref + E2,
                               col = dir_deriv), alpha = 0.7,
                           size=1) +  scale_colour_gradient2(low = muted("blue"),
                                                             mid = "white",
                                                             high = muted("red"),
                                                             midpoint = 0)



refs_sp2 <- refs %>% dplyr::filter(sp == "s2")
dd1_sp2<- tibble(angle = rep(seq(0.8, 2*pi, length = num_arrows), nrow(refs_sp2))) %>% 
  mutate(E1_ref = rep(refs_sp2$E1_ref, each = num_arrows),
         E2_ref = rep(refs_sp2$E2_ref, each = num_arrows)) %>%  
  full_join(refs_sp2) %>%  
  mutate(E1 = cos(angle) * radius,
         E2 = sin(angle) * radius,
         dir_deriv = E1 * pd_E1 + E2 * pd_E2, unit_vec_mag = sqrt(E1^2 + E2^2))

sp2 <- sp2 +  geom_segment(data = dd1_sp2,
                   aes(x = E1_ref, y = E2_ref,
                       xend = E1_ref + E1, yend = E2_ref + E2,
                       col = dir_deriv), alpha = 0.7,
                   size=1) +  scale_colour_gradient2(low = muted("blue"),
                                                     mid = "white",
                                                     high = muted("red"),
                                                     midpoint = 0)




refs_sp3 <- refs %>% dplyr::filter(sp == "s3")
dd1_sp3<- tibble(angle = rep(seq(0.8, 2*pi, length = num_arrows), nrow(refs_sp3))) %>% 
  mutate(E1_ref = rep(refs_sp3$E1_ref, each = num_arrows),
         E2_ref = rep(refs_sp3$E2_ref, each = num_arrows)) %>%  
  full_join(refs_sp3) %>%  
  mutate(E1 = cos(angle) * radius,
         E2 = sin(angle) * radius,
         dir_deriv = E1 * pd_E1 + E2 * pd_E2, unit_vec_mag = sqrt(E1^2 + E2^2))

sp3 <- sp3 +  geom_segment(data = dd1_sp3,
                    aes(x = E1_ref, y = E2_ref,
                        xend = E1_ref + E1, yend = E2_ref + E2,
                        col = dir_deriv), alpha = 0.7,
                    size=1) +  scale_colour_gradient2(low = muted("blue"),
                                                      mid = "white",
                                                      high = muted("red"),
                                                      midpoint = 0)



refs_sp4 <- refs %>% dplyr::filter(sp == "s4")
dd1_sp4<- tibble(angle = rep(seq(0.8, 2*pi, length = num_arrows), nrow(refs_sp4))) %>% 
  mutate(E1_ref = rep(refs_sp4$E1_ref, each = num_arrows),
         E2_ref = rep(refs_sp4$E2_ref, each = num_arrows)) %>%  
  full_join(refs_sp4) %>%  
  mutate(E1 = cos(angle) * radius,
         E2 = sin(angle) * radius,
         dir_deriv = E1 * pd_E1 + E2 * pd_E2, unit_vec_mag = sqrt(E1^2 + E2^2))

sp4 <- sp4 +  geom_segment(data = dd1_sp4,
                           aes(x = E1_ref, y = E2_ref,
                               xend = E1_ref + E1, yend = E2_ref + E2,
                               col = dir_deriv), alpha = 0.7,
                           size=1) +  scale_colour_gradient2(low = muted("blue"),
                                                             mid = "white",
                                                             high = muted("red"),
                                                             midpoint = 0)

ggarrange(sp1, sp2, sp3, sp4, ncol=2, nrow=2, common.legend = TRUE, legend="left")


### Mean across the surface

dd1_sp1 <- dd1_sp1 %>%  dplyr::mutate(mean_dd = mean(dir_deriv))

dd1_sp1$E1_ref_try <- 300
dd1_sp1$E2_ref_try <- 20




aver <- dd1_sp1[1,]
sp1 +  geom_segment(data = aver,
                    aes(x = E1_ref_try, y = E2_ref_try,
                        xend = E1_ref_try + E1, yend = E2_ref_try + E2,
                        col = "orange"), alpha = 0.7,
                    size=1) +  scale_colour_gradient2(low = muted("blue"),
                                                      mid = "white",
                                                      high = muted("red"),
                                                      midpoint = 0)









rmax<-max(rdiv_sim_nocor$rdiv); rmin<-1
dmax<-max(rdiv_sim_nocor$sign); dmin<-0
# community 1 


dd_plot <- ggplot(data = red_spp, mapping = aes(x = time, y = (dir_deriv), col = sp)) +
  theme_bw(base_size = 12)+
  geom_line()+
  labs(x = "time",y = "Directional derivative", tag = "a")+
  geom_hline(yintercept = 0, linetype= "dashed")

rdiv_plot <- ggplot(data = rdiv_sim_nocor, mapping = aes(x = time, y = rdiv)) +
  geom_line()+
  labs(x = "time",y = "Dissimilarity (derivatives)", tag = "b") +
  theme_bw(base_size = 12)+
  geom_hline(yintercept = mean(rdiv_sim_nocor$rdiv),lty=2) +
  geom_richtext(x = 7,
                mapping = aes(y = rmax - 0.005),
                size=4.5,
                label.color = NA,
                label = paste("RDiv^Mean =",paste0(round(mean(rdiv_sim_nocor$rdiv),digits = 2))))+
  lims(y = c(rmin,rmax))


sing_plot <- ggplot(data = rdiv_sim_nocor, mapping = aes(x = time, y = sign)) +
  geom_line()+
  labs(x = "time",y = "Divergence", tag = "c") +
  theme_bw(base_size = 12)+
  geom_hline(yintercept = mean(rdiv_sim_nocor$sign),lty=2) +
  geom_richtext(x = 7,
                mapping = aes(y = dmax - 0.1),
                size=4.5,
                label.color = NA,
                label = paste("RDiv^Mean =",paste0(round(mean(rdiv_sim_nocor$sign),digits = 2))))+
  lims(y = c(dmin,dmax))




(dd_plot ) / (rdiv_plot) / (sing_plot)





p1 <- p1E1 + xlab("Temperature (k)") + ylab("Salinity (ppt)") + labs(title = "te(Temperature, Salinity)") + theme_classic(base_size = 15)
p2 <- p2   + xlab("Temperature (k)") + labs(title = "te(Temperature, Salinity)")+ theme_classic(base_size = 15)
p3 <- p3   + xlab("Temperature (k)") + labs(title = "te(Temperature, Salinity)")+ theme_classic(base_size = 15)
p4 <- p4 + xlab("Temperature (k)") + ylab("Salinity (ppt)") + labs(title = "te(Temperature, Salinity)")+ theme_classic(base_size = 15)
p5 <- p5   + xlab("Salinity (ppt)") + labs(title = "te(Temperature, Salinity)")+ theme_classic(base_size = 15)
p6 <- p6   + xlab("Salinity (ppt)") + labs(title = "te(Temperature, Salinity)")+ theme_classic(base_size = 15)
