
## Low fluctuation in environmental change - negative correlation 
Increasing variance in z1 - negative correlation between env variables (E1 and E2) with *low* fluctuations

We repeat the same steps as in the previous section, but this time E1 and E2 have a very small variance.

(ref:plot2env) Environmental variables with *low* variation and negative correlation. (a) E1 changing over time. (b) E2 changing over time.
```{r plot2env, fig.cap='(ref:plot2env)', fig.align="center", fig.height=6, fig.width=12}

library(MASS)
sample_size <- 15                                       
sample_meanvector <- c(mean(E1_series),  mean(E2_series))                                   
sample_covariance_matrix <- matrix(c(1, -1,
                                     -1, 1),
                                   ncol = 2)

# create bivariate normal distribution
set.seed(65354)
refs <- mvrnorm(n = sample_size,
                mu = sample_meanvector, 
                Sigma = sample_covariance_matrix) %>%
  as_tibble() %>% 
  dplyr::rename(E1 = "V1", E2 = "V2") %>% 
  dplyr::mutate(time = seq.int(nrow(.)))

# print top of distribution
cor(refs)

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
```


```{r results='hide'}
# creating communities
# we are not going to plot spp responses, as they are the same as before.
s <- 4
# community 1 - low variation in z1. 
z1_range <- 1.5
z2_range <- 15

# create reproducible community
set.seed(2456)
comm <- get_community(s, z1_range, z2_range)

# Get directional dreivatives
dir_1 <- get_directionals(comm, new_data, refs)

# community 2 - intermediate variation in z1. 
z1_range <- 25
z2_range <- 15
# b1_range <- 0.009
# b2_range <- 0.001

# create reproducible community
set.seed(2456)
comm <- get_community(s, z1_range, z2_range)

# Get directional dreivatives
dir_2 <- get_directionals(comm, new_data, refs)

# community 3 - large variation in z1. 
z1_range <- 50
z2_range <- 15
# b1_range <- 0.02
# b2_range <- 0.001

# create reproducible community
set.seed(2456)
comm <- get_community(s, z1_range, z2_range)

# Get directional dreivatives
dir_3 <- get_directionals(comm, new_data, refs)


```


### Response diversity calculation 

(ref:RDplot2) Directional derivatives and response diversity for the three different communities. 
a, b, c. Species directional derivatives over time. d, e, f. Response diversity measured as similarity-based diversity metric. g, h, i. Response diversity measured as divergence (sign sensitive). 

```{r RDplot2, fig.cap='(ref:RDplot1)', fig.align="center", fig.height=10, fig.width=16}
# community 1  - calculation RD

# from long to wide
rdiv_1 <- dir_1 %>%
  spread( sp, dir_deriv)
rdiv_1[is.na(rdiv_1)] <- 0

# actual calculation for only the same species used above
rdiv_1$rdiv<-apply(rdiv_1[,c(4:7)], 1, resp_div, sign_sens = F)
rdiv_1$sign<-apply(rdiv_1[,c(4:7)], 1, resp_div, sign_sens = T)
rdiv_1$Med<-mean(rdiv_1$rdiv)
rdiv_1$Med_sing<-mean(rdiv_1$sign)


# community 2  - calculation RD
# from long to wide
rdiv_2 <- dir_2 %>%
  spread( sp, dir_deriv)
rdiv_2[is.na(rdiv_2)] <- 0

# actual calculation for only the same species used above
rdiv_2$rdiv<-apply(rdiv_2[,c(4:7)], 1, resp_div, sign_sens = F)
rdiv_2$sign<-apply(rdiv_2[,c(4:7)], 1, resp_div, sign_sens = T)
rdiv_2$Med<-mean(rdiv_2$rdiv)
rdiv_2$Med_sing<-mean(rdiv_2$sign)

# community 3  - calculation RD
# from long to wide
rdiv_3 <- dir_3 %>%
  spread( sp, dir_deriv)
rdiv_3[is.na(rdiv_3)] <- 0

# actual calculation for only the same species used above
rdiv_3$rdiv<-apply(rdiv_3[,c(4:7)], 1, resp_div, sign_sens = F)
rdiv_3$sign<-apply(rdiv_3[,c(4:7)], 1, resp_div, sign_sens = T)
rdiv_3$Med<-mean(rdiv_3$rdiv)
rdiv_3$Med_sing<-mean(rdiv_3$sign)


plot_RD(dir_1, dir_2, dir_3, rdiv_1, rdiv_2, rdiv_3)
```


## High fluctuation in environmental change - positive correlation
Increasing variance in z1 - positive correlation between env variables (E1 and E2) with *high* fluctuations

Same steps has before (3 communities with increasing diversity ib z1 while z2 is fixed), but the correlation between E1 and E1 is negative.


(ref:plot3env) Environmental variables with high variation and positive correlation. (a) E1 changing over time. (b) E2 changing over time.
```{r plot3env, fig.cap='(ref:plot3env)', fig.align="center", fig.height=6, fig.width=12}
library(MASS)
sample_size <- 15                                       
sample_meanvector <- c(mean(E1_series),  mean(E2_series))                                   
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
```




```{r}
s <- 4
# community 1 - low variation in z1. 
z1_range <- 1.5
z2_range <- 15

# create reproducible community
set.seed(2456)
comm <- get_community(s, z1_range, z2_range)

# Get directional dreivatives
dir_1 <- get_directionals(comm, new_data, refs)

# community 2 - intermediate variation in z1. 
z1_range <- 25
z2_range <- 15
# b1_range <- 0.009
# b2_range <- 0.001

# create reproducible community
set.seed(2456)
comm <- get_community(s, z1_range, z2_range)

# Get directional dreivatives
dir_2 <- get_directionals(comm, new_data, refs)


# community 3 - large variation in z1. 
z1_range <- 50
z2_range <- 15
# b1_range <- 0.02
# b2_range <- 0.001

# create reproducible community
set.seed(2456)
comm <- get_community(s, z1_range, z2_range)

# Get directional dreivatives
dir_3 <- get_directionals(comm, new_data, refs)


```


### Response diversity calculation 

(ref:RDplot3) Directional derivatives and response diversity for the three different communities. 
a, b, c. Species directional derivatives over time. d, e, f. Response diversity measured as similarity-based diversity metric. g, h, i. Response diversity measured as divergence (sign sensitive). 
```{r RDplot3, fig.cap='(ref:RDplot3)', fig.align="center", fig.height=10, fig.width=16}
# community 1  - calculation RD

# from long to wide
rdiv_1 <- dir_1 %>%
  spread( sp, dir_deriv)
rdiv_1[is.na(rdiv_1)] <- 0

# actual calculation for only the same species used above
rdiv_1$rdiv<-apply(rdiv_1[,c(4:7)], 1, resp_div, sign_sens = F)
rdiv_1$sign<-apply(rdiv_1[,c(4:7)], 1, resp_div, sign_sens = T)
rdiv_1$Med<-mean(rdiv_1$rdiv)
rdiv_1$Med_sing<-mean(rdiv_1$sign)


# community 2  - calculation RD
# from long to wide
rdiv_2 <- dir_2 %>%
  spread( sp, dir_deriv)
rdiv_2[is.na(rdiv_2)] <- 0

# actual calculation for only the same species used above
rdiv_2$rdiv<-apply(rdiv_2[,c(4:7)], 1, resp_div, sign_sens = F)
rdiv_2$sign<-apply(rdiv_2[,c(4:7)], 1, resp_div, sign_sens = T)
rdiv_2$Med<-mean(rdiv_2$rdiv)
rdiv_2$Med_sing<-mean(rdiv_2$sign)

# community 3  - calculation RD
# from long to wide
rdiv_3 <- dir_3 %>%
  spread( sp, dir_deriv)
rdiv_3[is.na(rdiv_3)] <- 0

# actual calculation for only the same species used above
rdiv_3$rdiv<-apply(rdiv_3[,c(4:7)], 1, resp_div, sign_sens = F)
rdiv_3$sign<-apply(rdiv_3[,c(4:7)], 1, resp_div, sign_sens = T)
rdiv_3$Med<-mean(rdiv_3$rdiv)
rdiv_3$Med_sing<-mean(rdiv_3$sign)



plot_RD(dir_1, dir_2, dir_3, rdiv_1, rdiv_2, rdiv_3)
```

## Low fluctuation in environmental change - positive correlation
Increasing variance in z1 - positive correlation between env variables (E1 and E2) with *low* fluctuations

We repeat the same steps as in the previous section, but this time E1 and E2 have a very small variance.

(ref:plot2env) Environmental variables with *low* variation and negative correlation. (a) E1 changing over time. (b) E2 changing over time.
```{r plot4env, fig.cap='(ref:plot4env)', fig.align="center", fig.height=6, fig.width=12}
library(MASS)
sample_size <- 15                                       
sample_meanvector <- c(mean(E1_series),  mean(E2_series))                                   
sample_covariance_matrix <- matrix(c(1, 1,
                                     1, 1),
                                   ncol = 2)

# create bivariate normal distribution
set.seed(5655)
refs <- mvrnorm(n = sample_size,
                mu = sample_meanvector, 
                Sigma = sample_covariance_matrix) %>%
  as_tibble() %>% 
  dplyr::rename(E1 = "V1", E2 = "V2") %>% 
  dplyr::mutate(time = seq.int(nrow(.)))

# print top of distribution
cor(refs)

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
```

```{r}
s <- 4
# community 1 - low variation in z1. 
z1_range <- 1.5
z2_range <- 15

# create reproducible community
set.seed(2456)
comm <- get_community(s, z1_range, z2_range)

# Get directional dreivatives
dir_1 <- get_directionals(comm, new_data, refs)

# community 2 - intermediate variation in z1. 
z1_range <- 25
z2_range <- 15
# b1_range <- 0.009
# b2_range <- 0.001

# create reproducible community
set.seed(2456)
comm <- get_community(s, z1_range, z2_range)

# Get directional dreivatives
dir_2 <- get_directionals(comm, new_data, refs)

# community 3 - large variation in z1. 
z1_range <- 50
z2_range <- 15
# b1_range <- 0.02
# b2_range <- 0.001

# create reproducible community
set.seed(2456)
comm <- get_community(s, z1_range, z2_range)

# Get directional dreivatives
dir_3 <- get_directionals(comm, new_data, refs)

```


### Response diversity calculation 

(ref:RDplot4) Directional derivatives and response diversity for the three different communities. 
a, b, c. Species directional derivatives over time. d, e, f. Response diversity measured as similarity-based diversity metric. g, h, i. Response diversity measured as divergence (sign sensitive). 
```{r RDplot4, fig.cap='(ref:RDplot4)', fig.align="center", fig.height=10, fig.width=16}

# community 1  - calculation RD

# from long to wide
rdiv_1 <- dir_1 %>%
  spread( sp, dir_deriv)
rdiv_1[is.na(rdiv_1)] <- 0

# actual calculation for only the same species used above
rdiv_1$rdiv<-apply(rdiv_1[,c(4:7)], 1, resp_div, sign_sens = F)
rdiv_1$sign<-apply(rdiv_1[,c(4:7)], 1, resp_div, sign_sens = T)
rdiv_1$Med<-mean(rdiv_1$rdiv)
rdiv_1$Med_sing<-mean(rdiv_1$sign)


# community 2  - calculation RD
# from long to wide
rdiv_2 <- dir_2 %>%
  spread( sp, dir_deriv)
rdiv_2[is.na(rdiv_2)] <- 0

# actual calculation for only the same species used above
rdiv_2$rdiv<-apply(rdiv_2[,c(4:7)], 1, resp_div, sign_sens = F)
rdiv_2$sign<-apply(rdiv_2[,c(4:7)], 1, resp_div, sign_sens = T)
rdiv_2$Med<-mean(rdiv_2$rdiv)
rdiv_2$Med_sing<-mean(rdiv_2$sign)

# community 3  - calculation RD
# from long to wide
rdiv_3 <- dir_3 %>%
  spread( sp, dir_deriv)
rdiv_3[is.na(rdiv_3)] <- 0

# actual calculation for only the same species used above
rdiv_3$rdiv<-apply(rdiv_3[,c(4:7)], 1, resp_div, sign_sens = F)
rdiv_3$sign<-apply(rdiv_3[,c(4:7)], 1, resp_div, sign_sens = T)
rdiv_3$Med<-mean(rdiv_3$rdiv)
rdiv_3$Med_sing<-mean(rdiv_3$sign)


plot_RD(dir_1, dir_2, dir_3, rdiv_1, rdiv_2, rdiv_3)
```


## High fluctuation in environmental change - no correlation
Increasing variance in z1 - no correlation between env variables (E1 and E2) with *high* fluctuations

Same steps has before (3 communities with increasing diversity ib z1 while z2 is fixed), but there is no correlation between E1 and E1


(ref:plot5env) Environmental variables with high variation and no correlation. (a) E1 changing over time. (b) E2 changing over time.
```{r plot5env, fig.cap='(ref:plot5env)', fig.align="center", fig.height=6, fig.width=12}
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
```

```{r}
# creating the communities
s <- 4
# community 1 - low variation in z1. 
z1_range <- 1.5
z2_range <- 15

# create reproducible community
set.seed(2456)
comm <- get_community(s, z1_range, z2_range)

# Get directional dreivatives
dir_1 <- get_directionals(comm, new_data, refs)

# community 2 - intermediate variation in z1. 
z1_range <- 25
z2_range <- 15
# b1_range <- 0.009
# b2_range <- 0.001

# create reproducible community
set.seed(2456)
comm <- get_community(s, z1_range, z2_range)

# Get directional dreivatives
dir_2 <- get_directionals(comm, new_data, refs)

# community 3 - large variation in z1. 
z1_range <- 50
z2_range <- 15
# b1_range <- 0.02
# b2_range <- 0.001

# create reproducible community
set.seed(2456)
comm <- get_community(s, z1_range, z2_range)

# Get directional dreivatives
dir_3 <- get_directionals(comm, new_data, refs)

```


### Response diversity calculation 

(ref:RDplot5) Directional derivatives and response diversity for the three different communities. 
a, b, c. Species directional derivatives over time. d, e, f. Response diversity measured as similarity-based diversity metric. g, h, i. Response diversity measured as divergence (sign sensitive). 
```{r RDplot5, fig.cap='(ref:RDplot5)', fig.align="center", fig.height=10, fig.width=16}
# community 1  - calculation RD

# from long to wide
rdiv_1 <- dir_1 %>%
  spread( sp, dir_deriv)
rdiv_1[is.na(rdiv_1)] <- 0

# actual calculation for only the same species used above
rdiv_1$rdiv<-apply(rdiv_1[,c(4:7)], 1, resp_div, sign_sens = F)
rdiv_1$sign<-apply(rdiv_1[,c(4:7)], 1, resp_div, sign_sens = T)
rdiv_1$Med<-mean(rdiv_1$rdiv)
rdiv_1$Med_sing<-mean(rdiv_1$sign)

# community 2  - calculation RD
# from long to wide
rdiv_2 <- dir_2 %>%
  spread( sp, dir_deriv)
rdiv_2[is.na(rdiv_2)] <- 0

# actual calculation for only the same species used above
rdiv_2$rdiv<-apply(rdiv_2[,c(4:7)], 1, resp_div, sign_sens = F)
rdiv_2$sign<-apply(rdiv_2[,c(4:7)], 1, resp_div, sign_sens = T)
rdiv_2$Med<-mean(rdiv_2$rdiv)
rdiv_2$Med_sing<-mean(rdiv_2$sign)

# community 3  - calculation RD
# from long to wide
rdiv_3 <- dir_3 %>%
  spread( sp, dir_deriv)
rdiv_3[is.na(rdiv_3)] <- 0

# actual calculation for only the same species used above
rdiv_3$rdiv<-apply(rdiv_3[,c(4:7)], 1, resp_div, sign_sens = F)
rdiv_3$sign<-apply(rdiv_3[,c(4:7)], 1, resp_div, sign_sens = T)
rdiv_3$Med<-mean(rdiv_3$rdiv)
rdiv_3$Med_sing<-mean(rdiv_3$sign)

plot_RD(dir_1, dir_2, dir_3, rdiv_1, rdiv_2, rdiv_3)
```


## Low fluctuation in environmental change - no correlation
Increasing variance in z1 - no correlation between env variables (E1 and E2) with *low* fluctuations

We repeat the same steps as in the previous section, but this time E1 and E2 have a very small variance.

(ref:plot6env) Environmental variables with *low* variation and negative correlation. (a) E1 changing over time. (b) E2 changing over time.
```{r plot6env, fig.cap='(ref:plot6env)', fig.align="center", fig.height=6, fig.width=12}

library(MASS)
sample_size <- 15                                       
sample_meanvector <- c(mean(E1_series),  mean(E2_series))                                   
sample_covariance_matrix <- matrix(c(1, 0.0,
                                     0.0, 1),
                                   ncol = 2)

# create bivariate normal distribution
set.seed(5655)
refs <- mvrnorm(n = sample_size,
                mu = sample_meanvector, 
                Sigma = sample_covariance_matrix) %>%
  as_tibble() %>% 
  dplyr::rename(E1 = "V1", E2 = "V2") %>% 
  dplyr::mutate(time = seq.int(nrow(.)))

# print top of distribution
cor(refs)

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
```


```{r}
# create communities

s <- 4
# community 1 - low variation in z1. 
z1_range <- 1.5
z2_range <- 15

# create reproducible community
set.seed(2456)
comm <- get_community(s, z1_range, z2_range)

# Get directional dreivatives
dir_1 <- get_directionals(comm, new_data, refs)

head(dir_1)


# community 2 - intermediate variation in z1. 
z1_range <- 25
z2_range <- 15
# b1_range <- 0.009
# b2_range <- 0.001

# create reproducible community
set.seed(2456)
comm <- get_community(s, z1_range, z2_range)

# Get directional dreivatives
dir_2 <- get_directionals(comm, new_data, refs)

head(dir_2)


# community 3 - large variation in z1. 
z1_range <- 50
z2_range <- 15
# b1_range <- 0.02
# b2_range <- 0.001

# create reproducible community
set.seed(2456)
comm <- get_community(s, z1_range, z2_range)

# Get directional dreivatives
dir_3 <- get_directionals(comm, new_data, refs)

```


### Response diversity calculation 

(ref:RDplot6) Directional derivatives and response diversity for the three different communities. 
a, b, c. Species directional derivatives over time. d, e, f. Response diversity measured as similarity-based diversity metric. g, h, i. Response diversity measured as divergence (sign sensitive). 
```{r RDplot6, fig.cap='(ref:RDplot6)', fig.align="center", fig.height=10, fig.width=16}
# community 1  - calculation RD

# from long to wide
rdiv_1 <- dir_1 %>%
  spread( sp, dir_deriv)
rdiv_1[is.na(rdiv_1)] <- 0

# actual calculation for only the same species used above
rdiv_1$rdiv<-apply(rdiv_1[,c(4:7)], 1, resp_div, sign_sens = F)
rdiv_1$sign<-apply(rdiv_1[,c(4:7)], 1, resp_div, sign_sens = T)
rdiv_1$Med<-mean(rdiv_1$rdiv)
rdiv_1$Med_sing<-mean(rdiv_1$sign)

# community 2  - calculation RD
# from long to wide
rdiv_2 <- dir_2 %>%
  spread( sp, dir_deriv)
rdiv_2[is.na(rdiv_2)] <- 0

# actual calculation for only the same species used above
rdiv_2$rdiv<-apply(rdiv_2[,c(4:7)], 1, resp_div, sign_sens = F)
rdiv_2$sign<-apply(rdiv_2[,c(4:7)], 1, resp_div, sign_sens = T)
rdiv_2$Med<-mean(rdiv_2$rdiv)
rdiv_2$Med_sing<-mean(rdiv_2$sign)

# community 3  - calculation RD
# from long to wide
rdiv_3 <- dir_3 %>%
  spread( sp, dir_deriv)
rdiv_3[is.na(rdiv_3)] <- 0

# actual calculation for only the same species used above
rdiv_3$rdiv<-apply(rdiv_3[,c(4:7)], 1, resp_div, sign_sens = F)
rdiv_3$sign<-apply(rdiv_3[,c(4:7)], 1, resp_div, sign_sens = T)
rdiv_3$Med<-mean(rdiv_3$rdiv)
rdiv_3$Med_sing<-mean(rdiv_3$sign)

plot_RD(dir_1, dir_2, dir_3, rdiv_1, rdiv_2, rdiv_3)
```



# Variation in diversity in z1 and z2
Now we do the same, but having also z2 changing together with z1 with positive correlation. So, diversity in z1 and z2 are increasing gradually from low to high in the 3 communities.
I expect RD to increase when both z1 and z2 are high and there is negative correlation in env change (E1 increases when E2 decreases). 


## High fluctuation in environmental change - negative correlation
Increasing variance in z1 and z2 - negative correlation between env variables (E1 and E2) with *high* fluctuations
We now create 3 communities. Community 1 is characterised by low diversity in z1 and z2, community 2 has medium diversity in z1 and z2, and community 3 has high diversity in z1 and z2. 

E1 and E2 have high fluctuations and negative correlation. 

(ref:plot7env) Environmental variables with *high* variation and negative correlation. (a) E1 changing over time. (b) E2 changing over time.
```{r plot7env, fig.cap='(ref:plot7env)', fig.align="center", fig.height=6, fig.width=12}
# Specify the correlation E1 and E2 
#negative in this example and high variation in E1 and E2 

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

```

### Community 1 - low diversity in z1 and z2. 

(ref:plotcomm1.2_neg) Community 1, single species responses. (a) Species responses to the gradient of E1. Different colour lines show the dependency of the rate to the second environmental variable (E2). (b) Species responses to the gradient of E2. Different colour lines show the dependency of the rate to the second environmental variable (E1).
```{r plotcomm1.2_neg, fig.cap='(ref:plotcomm1.2_neg)', fig.align="center", fig.height=6, fig.width=12}
s <- 4
# community 1 - low variation in z1. 
z1_range <- 1.5
z2_range <- 1.5

# create reproducible community
set.seed(2456)
comm <- get_community(s, z1_range, z2_range)

# Visualise spp performances 
performance_comm1 <- multi_performance_plot(comm)
performance_comm1
# Get directional dreivatives
dir_1 <- get_directionals(comm, new_data, refs)

```


### Community 2 - medium diversity in z1 and z2. 

(ref:plotcomm2.1_neg) Community 2, single species responses. (a) Species responses to the gradient of E1. Different colour lines show the dependency of the rate to the second environmental variable (E2). (b) Species responses to the gradient of E2. Different colour lines show the dependency of the rate to the second environmental variable (E1).
```{r plotcomm2.1_neg, fig.cap='(ref:plotcomm2.1_neg)', fig.align="center", fig.height=6, fig.width=12}
z1_range <- 25
z2_range <- 25
# b1_range <- 0.009
# b2_range <- 0.001

# create reproducible community
set.seed(2456)
comm <- get_community(s, z1_range, z2_range)

# Visualise spp performances 
performance_comm2 <- multi_performance_plot(comm)
performance_comm2
# Get directional dreivatives
dir_2 <- get_directionals(comm, new_data, refs)

```



### Community 3 - high diversity in z1 and z2. 

(ref:plotcomm3.2_neg) Community 3, single species responses. (a) Species responses to the gradient of E1. Different colour lines show the dependency of the rate to the second environmental variable (E2). (b) Species responses to the gradient of E2. Different colour lines show the dependency of the rate to the second environmental variable (E1).
```{r plotcomm3.2_neg, fig.cap='(ref:plotcomm3.2_neg)', fig.align="center", fig.height=6, fig.width=12}
z1_range <- 50
z2_range <- 50
# b1_range <- 0.009
# b2_range <- 0.001

# create reproducible community
set.seed(2456)
comm <- get_community(s, z1_range, z2_range)

# Visualise spp performances 
performance_comm3 <- multi_performance_plot(comm)
performance_comm3
# Get directional dreivatives
dir_3 <- get_directionals(comm, new_data, refs)
```


### Response diversity calculation 

(ref:RDplot7) Directional derivatives and response diversity for the three different communities. 
a, b, c. Species directional derivatives over time. d, e, f. Response diversity measured as similarity-based diversity metric. g, h, i. Response diversity measured as divergence (sign sensitive). 
```{r RDplot7, fig.cap='(ref:RDplot7)', fig.align="center", fig.height=10, fig.width=16}

# community 1  - calculation RD

# from long to wide
rdiv_1 <- dir_1 %>%
  spread( sp, dir_deriv)
rdiv_1[is.na(rdiv_1)] <- 0

# actual calculation for only the same species used above
rdiv_1$rdiv<-apply(rdiv_1[,c(4:7)], 1, resp_div, sign_sens = F)
rdiv_1$sign<-apply(rdiv_1[,c(4:7)], 1, resp_div, sign_sens = T)
rdiv_1$Med<-mean(rdiv_1$rdiv)
rdiv_1$Med_sing<-mean(rdiv_1$sign)

# community 2  - calculation RD
# from long to wide
rdiv_2 <- dir_2 %>%
  spread( sp, dir_deriv)
rdiv_2[is.na(rdiv_2)] <- 0

# actual calculation for only the same species used above
rdiv_2$rdiv<-apply(rdiv_2[,c(4:7)], 1, resp_div, sign_sens = F)
rdiv_2$sign<-apply(rdiv_2[,c(4:7)], 1, resp_div, sign_sens = T)
rdiv_2$Med<-mean(rdiv_2$rdiv)
rdiv_2$Med_sing<-mean(rdiv_2$sign)

# community 3  - calculation RD
# from long to wide
rdiv_3 <- dir_3 %>%
  spread( sp, dir_deriv)
rdiv_3[is.na(rdiv_3)] <- 0

# actual calculation for only the same species used above
rdiv_3$rdiv<-apply(rdiv_3[,c(4:7)], 1, resp_div, sign_sens = F)
rdiv_3$sign<-apply(rdiv_3[,c(4:7)], 1, resp_div, sign_sens = T)
rdiv_3$Med<-mean(rdiv_3$rdiv)
rdiv_3$Med_sing<-mean(rdiv_3$sign)

plot_RD(dir_1, dir_2, dir_3, rdiv_1, rdiv_2, rdiv_3)
```


## Low fluctuation in environmental change - negative correlation 
Increasing variance in z1 and z2 - negative correlation between env variables (E1 and E2) with *low* fluctuations

We repeat the same steps as in the previous section, but this time E1 and E2 have a very small variance.

(ref:plot8env) Environmental variables with *low* variation and negative correlation. (a) E1 changing over time. (b) E2 changing over time.
```{r plot8env, fig.cap='(ref:plot8env)', fig.align="center", fig.height=6, fig.width=12}

library(MASS)
sample_size <- 15                                       
sample_meanvector <- c(mean(E1_series),  mean(E2_series))                                   
sample_covariance_matrix <- matrix(c(1, -1,
                                     -1, 1),
                                   ncol = 2)

# create bivariate normal distribution
set.seed(65354)
refs <- mvrnorm(n = sample_size,
                mu = sample_meanvector, 
                Sigma = sample_covariance_matrix) %>%
  as_tibble() %>% 
  dplyr::rename(E1 = "V1", E2 = "V2") %>% 
  dplyr::mutate(time = seq.int(nrow(.)))

# print top of distribution
cor(refs)

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
```

```{r}
# create communities
s <- 4
# community 1 - low variation in z1. 
z1_range <- 1.5
z2_range <- 1.5

# create reproducible community
set.seed(2456)
comm <- get_community(s, z1_range, z2_range)

# Get directional dreivatives
dir_1 <- get_directionals(comm, new_data, refs)

# community 2 - intermediate variation in z1. 
z1_range <- 25
z2_range <- 25
# b1_range <- 0.009
# b2_range <- 0.001

# create reproducible community
set.seed(2456)
comm <- get_community(s, z1_range, z2_range)

# Get directional dreivatives
dir_2 <- get_directionals(comm, new_data, refs)

# community 3 - large variation in z1. 
z1_range <- 50
z2_range <- 50
# b1_range <- 0.02
# b2_range <- 0.001

# create reproducible community
set.seed(2456)
comm <- get_community(s, z1_range, z2_range)

# Get directional dreivatives
dir_3 <- get_directionals(comm, new_data, refs)

```


### Response diversity calculation 

(ref:RDplot8) Directional derivatives and response diversity for the three different communities. 
a, b, c. Species directional derivatives over time. d, e, f. Response diversity measured as similarity-based diversity metric. g, h, i. Response diversity measured as divergence (sign sensitive). 
```{r RDplot8, fig.cap='(ref:RDplot8)', fig.align="center", fig.height=10, fig.width=16}
# community 1  - calculation RD

# from long to wide
rdiv_1 <- dir_1 %>%
  spread( sp, dir_deriv)
rdiv_1[is.na(rdiv_1)] <- 0

# actual calculation for only the same species used above
rdiv_1$rdiv<-apply(rdiv_1[,c(4:7)], 1, resp_div, sign_sens = F)
rdiv_1$sign<-apply(rdiv_1[,c(4:7)], 1, resp_div, sign_sens = T)
rdiv_1$Med<-mean(rdiv_1$rdiv)
rdiv_1$Med_sing<-mean(rdiv_1$sign)

# community 2  - calculation RD
# from long to wide
rdiv_2 <- dir_2 %>%
  spread( sp, dir_deriv)
rdiv_2[is.na(rdiv_2)] <- 0

# actual calculation for only the same species used above
rdiv_2$rdiv<-apply(rdiv_2[,c(4:7)], 1, resp_div, sign_sens = F)
rdiv_2$sign<-apply(rdiv_2[,c(4:7)], 1, resp_div, sign_sens = T)
rdiv_2$Med<-mean(rdiv_2$rdiv)
rdiv_2$Med_sing<-mean(rdiv_2$sign)

# community 3  - calculation RD
# from long to wide
rdiv_3 <- dir_3 %>%
  spread( sp, dir_deriv)
rdiv_3[is.na(rdiv_3)] <- 0

# actual calculation for only the same species used above
rdiv_3$rdiv<-apply(rdiv_3[,c(4:7)], 1, resp_div, sign_sens = F)
rdiv_3$sign<-apply(rdiv_3[,c(4:7)], 1, resp_div, sign_sens = T)
rdiv_3$Med<-mean(rdiv_3$rdiv)
rdiv_3$Med_sing<-mean(rdiv_3$sign)

plot_RD(dir_1, dir_2, dir_3, rdiv_1, rdiv_2, rdiv_3)
```



## High fluctuation in environmental change - positive correlation
Increasing variance in z1 and z2 - positive correlation between env variables (E1 and E2) with *high* fluctuations

Same steps has before (3 communities with increasing diversity in z1 and z2), but there is positive correlation between E1 and E1


(ref:plot9env) Environmental variables with high variation and positive correlation. (a) E1 changing over time. (b) E2 changing over time.
```{r plot9env, fig.cap='(ref:plot9env)', fig.align="center", fig.height=6, fig.width=12}
library(MASS)
sample_size <- 15                                       
sample_meanvector <- c(mean(E1_series),  mean(E2_series))                                   
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

# print top of distribution
cor(refs)

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
```



```{r}

# create communities
s <- 4
# community 1 - low variation in z1. 
z1_range <- 1.5
z2_range <- 1.5

# create reproducible community
set.seed(2456)
comm <- get_community(s, z1_range, z2_range)

# Get directional dreivatives
dir_1 <- get_directionals(comm, new_data, refs)

# community 2 - intermediate variation in z1. 
z1_range <- 25
z2_range <- 25
# b1_range <- 0.009
# b2_range <- 0.001

# create reproducible community
set.seed(2456)
comm <- get_community(s, z1_range, z2_range)

# Get directional dreivatives
dir_2 <- get_directionals(comm, new_data, refs)

# community 3 - large variation in z1. 
z1_range <- 50
z2_range <- 50
# b1_range <- 0.02
# b2_range <- 0.001

# create reproducible community
set.seed(2456)
comm <- get_community(s, z1_range, z2_range)

# Get directional dreivatives
dir_3 <- get_directionals(comm, new_data, refs)

```




### Response diversity calculation 

(ref:RDplot9) Directional derivatives and response diversity for the three different communities. 
a, b, c. Species directional derivatives over time. d, e, f. Response diversity measured as similarity-based diversity metric. g, h, i. Response diversity measured as divergence (sign sensitive). 
```{r RDplot9, fig.cap='(ref:RDplot9)', fig.align="center", fig.height=10, fig.width=16}

# community 1  - calculation RD

# from long to wide
rdiv_1 <- dir_1 %>%
  spread( sp, dir_deriv)
rdiv_1[is.na(rdiv_1)] <- 0

# actual calculation for only the same species used above
rdiv_1$rdiv<-apply(rdiv_1[,c(4:7)], 1, resp_div, sign_sens = F)
rdiv_1$sign<-apply(rdiv_1[,c(4:7)], 1, resp_div, sign_sens = T)
rdiv_1$Med<-mean(rdiv_1$rdiv)
rdiv_1$Med_sing<-mean(rdiv_1$sign)

# community 2  - calculation RD
# from long to wide
rdiv_2 <- dir_2 %>%
  spread( sp, dir_deriv)
rdiv_2[is.na(rdiv_2)] <- 0

# actual calculation for only the same species used above
rdiv_2$rdiv<-apply(rdiv_2[,c(4:7)], 1, resp_div, sign_sens = F)
rdiv_2$sign<-apply(rdiv_2[,c(4:7)], 1, resp_div, sign_sens = T)
rdiv_2$Med<-mean(rdiv_2$rdiv)
rdiv_2$Med_sing<-mean(rdiv_2$sign)

# community 3  - calculation RD
# from long to wide
rdiv_3 <- dir_3 %>%
  spread( sp, dir_deriv)
rdiv_3[is.na(rdiv_3)] <- 0

# actual calculation for only the same species used above
rdiv_3$rdiv<-apply(rdiv_3[,c(4:7)], 1, resp_div, sign_sens = F)
rdiv_3$sign<-apply(rdiv_3[,c(4:7)], 1, resp_div, sign_sens = T)
rdiv_3$Med<-mean(rdiv_3$rdiv)
rdiv_3$Med_sing<-mean(rdiv_3$sign)

plot_RD(dir_1, dir_2, dir_3, rdiv_1, rdiv_2, rdiv_3)
```


## Low fluctuation in environmental change - positive correlation 
Increasing variance in z1 and z2 - positive correlation between env variables (E1 and E2) with *low* fluctuations

We repeat the same steps as in the previous section, but this time E1 and E2 have a very small variance.

(ref:plot10env) Environmental variables with *low* variation and positive correlation. (a) E1 changing over time. (b) E2 changing over time.
```{r plot10env, fig.cap='(ref:plot10env)', fig.align="center", fig.height=6, fig.width=12}

library(MASS)
sample_size <- 15                                       
sample_meanvector <- c(mean(E1_series),  mean(E2_series))                                   
sample_covariance_matrix <- matrix(c(1, 1,
                                     1, 1),
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
```


```{r}
s <- 4
# community 1 - low variation in z1. 
z1_range <- 1.5
z2_range <- 1.5

# create reproducible community
set.seed(2456)
comm <- get_community(s, z1_range, z2_range)

# Get directional dreivatives
dir_1 <- get_directionals(comm, new_data, refs)

head(dir_1)


# community 2 - intermediate variation in z1. 
z1_range <- 25
z2_range <- 25
# b1_range <- 0.009
# b2_range <- 0.001

# create reproducible community
set.seed(2456)
comm <- get_community(s, z1_range, z2_range)

# Get directional dreivatives
dir_2 <- get_directionals(comm, new_data, refs)

# community 3 - large variation in z1. 
z1_range <- 50
z2_range <- 50
# b1_range <- 0.02
# b2_range <- 0.001

# create reproducible community
set.seed(2456)
comm <- get_community(s, z1_range, z2_range)

# Get directional dreivatives
dir_3 <- get_directionals(comm, new_data, refs)
```



### Response diversity calculation 

(ref:RDplot10) Directional derivatives and response diversity for the three different communities. 
a, b, c. Species directional derivatives over time. d, e, f. Response diversity measured as similarity-based diversity metric. g, h, i. Response diversity measured as divergence (sign sensitive). 
```{r RDplot10, fig.cap='(ref:RDplot10)', fig.align="center", fig.height=10, fig.width=16}
# community 1  - calculation RD

# from long to wide
rdiv_1 <- dir_1 %>%
  spread( sp, dir_deriv)
rdiv_1[is.na(rdiv_1)] <- 0

# actual calculation for only the same species used above
rdiv_1$rdiv<-apply(rdiv_1[,c(4:7)], 1, resp_div, sign_sens = F)
rdiv_1$sign<-apply(rdiv_1[,c(4:7)], 1, resp_div, sign_sens = T)
rdiv_1$Med<-mean(rdiv_1$rdiv)
rdiv_1$Med_sing<-mean(rdiv_1$sign)

# community 2  - calculation RD
# from long to wide
rdiv_2 <- dir_2 %>%
  spread( sp, dir_deriv)
rdiv_2[is.na(rdiv_2)] <- 0

# actual calculation for only the same species used above
rdiv_2$rdiv<-apply(rdiv_2[,c(4:7)], 1, resp_div, sign_sens = F)
rdiv_2$sign<-apply(rdiv_2[,c(4:7)], 1, resp_div, sign_sens = T)
rdiv_2$Med<-mean(rdiv_2$rdiv)
rdiv_2$Med_sing<-mean(rdiv_2$sign)

# community 3  - calculation RD
# from long to wide
rdiv_3 <- dir_3 %>%
  spread( sp, dir_deriv)
rdiv_3[is.na(rdiv_3)] <- 0

# actual calculation for only the same species used above
rdiv_3$rdiv<-apply(rdiv_3[,c(4:7)], 1, resp_div, sign_sens = F)
rdiv_3$sign<-apply(rdiv_3[,c(4:7)], 1, resp_div, sign_sens = T)
rdiv_3$Med<-mean(rdiv_3$rdiv)
rdiv_3$Med_sing<-mean(rdiv_3$sign)

plot_RD(dir_1, dir_2, dir_3, rdiv_1, rdiv_2, rdiv_3)
```


## High fluctuation in environmental change - no correlation
Increasing variance in z1 and z2 - no correlation between env variables (E1 and E2) with *high* fluctuations

Same steps has before (3 communities with increasing diversity in z1 and z2), but there is no correlation between E1 and E1


(ref:plot11env) Environmental variables with high variation and no correlation. (a) E1 changing over time. (b) E2 changing over time.
```{r plot11env, fig.cap='(ref:plot11env)', fig.align="center", fig.height=6, fig.width=12}
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

# print top of distribution
cor(refs)

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
```


```{r}
# creating communities
s <- 4
# community 1 - low variation in z1. 
z1_range <- 1.5
z2_range <- 1.5

# create reproducible community
set.seed(2456)
comm <- get_community(s, z1_range, z2_range)

# Get directional dreivatives
dir_1 <- get_directionals(comm, new_data, refs)

head(dir_1)


# community 2 - intermediate variation in z1. 
z1_range <- 25
z2_range <- 25
# b1_range <- 0.009
# b2_range <- 0.001

# create reproducible community
set.seed(2456)
comm <- get_community(s, z1_range, z2_range)

# Get directional dreivatives
dir_2 <- get_directionals(comm, new_data, refs)

head(dir_2)


# community 3 - large variation in z1. 
z1_range <- 50
z2_range <- 50
# b1_range <- 0.02
# b2_range <- 0.001

# create reproducible community
set.seed(2456)
comm <- get_community(s, z1_range, z2_range)

# Get directional dreivatives
dir_3 <- get_directionals(comm, new_data, refs)

```


### Response diversity calculation 

(ref:RDplot11) Directional derivatives and response diversity for the three different communities. 
a, b, c. Species directional derivatives over time. d, e, f. Response diversity measured as similarity-based diversity metric. g, h, i. Response diversity measured as divergence (sign sensitive). 
```{r RDplot11, fig.cap='(ref:RDplot11)', fig.align="center", fig.height=10, fig.width=16}

# community 1  - calculation RD

# from long to wide
rdiv_1 <- dir_1 %>%
  spread( sp, dir_deriv)
rdiv_1[is.na(rdiv_1)] <- 0

# actual calculation for only the same species used above
rdiv_1$rdiv<-apply(rdiv_1[,c(4:7)], 1, resp_div, sign_sens = F)
rdiv_1$sign<-apply(rdiv_1[,c(4:7)], 1, resp_div, sign_sens = T)
rdiv_1$Med<-mean(rdiv_1$rdiv)
rdiv_1$Med_sing<-mean(rdiv_1$sign)

# community 2  - calculation RD
# from long to wide
rdiv_2 <- dir_2 %>%
  spread( sp, dir_deriv)
rdiv_2[is.na(rdiv_2)] <- 0

# actual calculation for only the same species used above
rdiv_2$rdiv<-apply(rdiv_2[,c(4:7)], 1, resp_div, sign_sens = F)
rdiv_2$sign<-apply(rdiv_2[,c(4:7)], 1, resp_div, sign_sens = T)
rdiv_2$Med<-mean(rdiv_2$rdiv)
rdiv_2$Med_sing<-mean(rdiv_2$sign)

# community 3  - calculation RD
# from long to wide
rdiv_3 <- dir_3 %>%
  spread( sp, dir_deriv)
rdiv_3[is.na(rdiv_3)] <- 0

# actual calculation for only the same species used above
rdiv_3$rdiv<-apply(rdiv_3[,c(4:7)], 1, resp_div, sign_sens = F)
rdiv_3$sign<-apply(rdiv_3[,c(4:7)], 1, resp_div, sign_sens = T)
rdiv_3$Med<-mean(rdiv_3$rdiv)
rdiv_3$Med_sing<-mean(rdiv_3$sign)

plot_RD(dir_1, dir_2, dir_3, rdiv_1, rdiv_2, rdiv_3)
```



## Low fluctuation in environmental change - no correlation
Increasing variance in z1 and z2 - no correlation between env variables (E1 and E2) with *low* fluctuations

We repeat the same steps as in the previous section, but this time E1 and E2 have a very small variance.

(ref:plot12env) Environmental variables with *low* variation and negative correlation. (a) E1 changing over time. (b) E2 changing over time.
```{r plot12env, fig.cap='(ref:plot12env)', fig.align="center", fig.height=6, fig.width=12}
library(MASS)
sample_size <- 15                                       
sample_meanvector <- c(mean(E1_series),  mean(E2_series))                                   
sample_covariance_matrix <- matrix(c(1, 0.0,
                                     0.0, 1),
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
```

```{r}

# create communities
s <- 4
# community 1 - low variation in z1. 
z1_range <- 1.5
z2_range <- 1.5

# create reproducible community
set.seed(2456)
comm <- get_community(s, z1_range, z2_range)

# Get directional dreivatives
dir_1 <- get_directionals(comm, new_data, refs)


# community 2 - intermediate variation in z1. 
z1_range <- 25
z2_range <- 25
# b1_range <- 0.009
# b2_range <- 0.001

# create reproducible community
set.seed(2456)
comm <- get_community(s, z1_range, z2_range)

# Get directional dreivatives
dir_2 <- get_directionals(comm, new_data, refs)

head(dir_2)


# community 3 - large variation in z1. 
z1_range <- 50
z2_range <- 50
# b1_range <- 0.02
# b2_range <- 0.001

# create reproducible community
set.seed(2456)
comm <- get_community(s, z1_range, z2_range)

# Get directional dreivatives
dir_3 <- get_directionals(comm, new_data, refs)


```

### Response diversity calculation 

(ref:RDplot12) Directional derivatives and response diversity for the three different communities. 
a, b, c. Species directional derivatives over time. d, e, f. Response diversity measured as similarity-based diversity metric. g, h, i. Response diversity measured as divergence (sign sensitive). 
```{r RDplot12, fig.cap='(ref:RDplot12)', fig.align="center", fig.height=10, fig.width=16}
# community 1  - calculation RD

# from long to wide
rdiv_1 <- dir_1 %>%
  spread( sp, dir_deriv)
rdiv_1[is.na(rdiv_1)] <- 0

# actual calculation for only the same species used above
rdiv_1$rdiv<-apply(rdiv_1[,c(4:7)], 1, resp_div, sign_sens = F)
rdiv_1$sign<-apply(rdiv_1[,c(4:7)], 1, resp_div, sign_sens = T)
rdiv_1$Med<-mean(rdiv_1$rdiv)
rdiv_1$Med_sing<-mean(rdiv_1$sign)

# community 2  - calculation RD
# from long to wide
rdiv_2 <- dir_2 %>%
  spread( sp, dir_deriv)
rdiv_2[is.na(rdiv_2)] <- 0

# actual calculation for only the same species used above
rdiv_2$rdiv<-apply(rdiv_2[,c(4:7)], 1, resp_div, sign_sens = F)
rdiv_2$sign<-apply(rdiv_2[,c(4:7)], 1, resp_div, sign_sens = T)
rdiv_2$Med<-mean(rdiv_2$rdiv)
rdiv_2$Med_sing<-mean(rdiv_2$sign)

# community 3  - calculation RD
# from long to wide
rdiv_3 <- dir_3 %>%
  spread( sp, dir_deriv)
rdiv_3[is.na(rdiv_3)] <- 0

# actual calculation for only the same species used above
rdiv_3$rdiv<-apply(rdiv_3[,c(4:7)], 1, resp_div, sign_sens = F)
rdiv_3$sign<-apply(rdiv_3[,c(4:7)], 1, resp_div, sign_sens = T)
rdiv_3$Med<-mean(rdiv_3$rdiv)
rdiv_3$Med_sing<-mean(rdiv_3$sign)


plot_RD(dir_1, dir_2, dir_3, rdiv_1, rdiv_2, rdiv_3)
```


# Variation in diversity in z1 and z2
Now we do the same, but having also z2 changing together with z1 with negative correlation. So, diversity in z1 is increasing gradually from low to high in the 3 communities, while z2 decreases gradually in the 3 communities.

I expect RD to be the highest when either z1 or z2 is high and there is negative correlation in env change (E1 increases when E2 decreases). 


## High fluctuation in environmental change - negative correlation
Increasing variance in z1 and z2 - negative correlation between env variables (E1 and E2) with *high* fluctuations
We now create 3 communities. Community 1 is characterised by low diversity in z1 and high diversity in z2, community 2 has medium diversity in z1 and z2, and community 3 has high diversity in z1 and low in z2. 

E1 and E2 have high fluctuations and negative correlation. 

(ref:plot13env) Environmental variables with *high* variation and negative correlation. (a) E1 changing over time. (b) E2 changing over time.
```{r plot13env, fig.cap='(ref:plot13env)', fig.align="center", fig.height=6, fig.width=12}
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
```

### Community 1 - low diversity in z1 and high diversity in z2. 

(ref:plotcomm1.3_neg) Community 1, single species responses. (a) Species responses to the gradient of E1. Different colour lines show the dependency of the rate to the second environmental variable (E2). (b) Species responses to the gradient of E2. Different colour lines show the dependency of the rate to the second environmental variable (E1).
```{r plotcomm1.3_neg, fig.cap='(ref:plotcomm1.3_neg)', fig.align="center", fig.height=6, fig.width=12}
s <- 4
# community 1 - low variation in z1 high in z2. 
z1_range <- 1.5
z2_range <- 50

# create reproducible community
set.seed(2456)
comm <- get_community(s, z1_range, z2_range)
# Visualise spp performances 
performance_comm1 <- multi_performance_plot(comm)
performance_comm1

# Get directional dreivatives
dir_1 <- get_directionals(comm, new_data, refs)
```


### Community 2 - medium diversity in z1 and z2. 

(ref:plotcomm2.3_neg) Community 2, single species responses. (a) Species responses to the gradient of E1. Different colour lines show the dependency of the rate to the second environmental variable (E2). (b) Species responses to the gradient of E2. Different colour lines show the dependency of the rate to the second environmental variable (E1).
```{r plotcomm2.3_neg, fig.cap='(ref:plotcomm2.3_neg)', fig.align="center", fig.height=6, fig.width=12}
z1_range <- 25
z2_range <- 25
# b1_range <- 0.009
# b2_range <- 0.001

# create reproducible community
set.seed(2456)
comm <- get_community(s, z1_range, z2_range)

# Visualise spp performances 
performance_comm2 <- multi_performance_plot(comm)
performance_comm2
# Get directional dreivatives
dir_2 <- get_directionals(comm, new_data, refs)
```

### Community 3 - high diversity in z1 and low in z2. 

(ref:plotcomm3.3_neg) Community 3, single species responses. (a) Species responses to the gradient of E1. Different colour lines show the dependency of the rate to the second environmental variable (E2). (b) Species responses to the gradient of E2. Different colour lines show the dependency of the rate to the second environmental variable (E1).
```{r plotcomm3.3_neg, fig.cap='(ref:plotcomm3.3_neg)', fig.align="center", fig.height=6, fig.width=12}
z1_range <- 50
z2_range <- 1.5
# b1_range <- 0.009
# b2_range <- 0.001

# create reproducible community
set.seed(2456)
comm <- get_community(s, z1_range, z2_range)

# Visualise spp performances 
performance_comm3 <- multi_performance_plot(comm)
performance_comm3
# Get directional dreivatives
dir_3 <- get_directionals(comm, new_data, refs)
```

### Response diversity calculation 

(ref:RDplot13) Directional derivatives and response diversity for the three different communities. 
a, b, c. Species directional derivatives over time. d, e, f. Response diversity measured as similarity-based diversity metric. g, h, i. Response diversity measured as divergence (sign sensitive). 
```{r RDplot13, fig.cap='(ref:RDplot13)', fig.align="center", fig.height=10, fig.width=16}
# community 1  - calculation RD

# from long to wide
rdiv_1 <- dir_1 %>%
  spread( sp, dir_deriv)
rdiv_1[is.na(rdiv_1)] <- 0

# actual calculation for only the same species used above
rdiv_1$rdiv<-apply(rdiv_1[,c(4:7)], 1, resp_div, sign_sens = F)
rdiv_1$sign<-apply(rdiv_1[,c(4:7)], 1, resp_div, sign_sens = T)
rdiv_1$Med<-mean(rdiv_1$rdiv)
rdiv_1$Med_sing<-mean(rdiv_1$sign)

# community 2  - calculation RD
# from long to wide
rdiv_2 <- dir_2 %>%
  spread( sp, dir_deriv)
rdiv_2[is.na(rdiv_2)] <- 0

# actual calculation for only the same species used above
rdiv_2$rdiv<-apply(rdiv_2[,c(4:7)], 1, resp_div, sign_sens = F)
rdiv_2$sign<-apply(rdiv_2[,c(4:7)], 1, resp_div, sign_sens = T)
rdiv_2$Med<-mean(rdiv_2$rdiv)
rdiv_2$Med_sing<-mean(rdiv_2$sign)

# community 3  - calculation RD
# from long to wide
rdiv_3 <- dir_3 %>%
  spread( sp, dir_deriv)
rdiv_3[is.na(rdiv_3)] <- 0

# actual calculation for only the same species used above
rdiv_3$rdiv<-apply(rdiv_3[,c(4:7)], 1, resp_div, sign_sens = F)
rdiv_3$sign<-apply(rdiv_3[,c(4:7)], 1, resp_div, sign_sens = T)
rdiv_3$Med<-mean(rdiv_3$rdiv)
rdiv_3$Med_sing<-mean(rdiv_3$sign)

plot_RD(dir_1, dir_2, dir_3, rdiv_1, rdiv_2, rdiv_3)
```



## Low fluctuation in environmental change - negative correlation 
Increasing variance in z1 and z2 - negative correlation between env variables (E1 and E2) with *low* fluctuations

We repeat the same steps as in the previous section, but this time E1 and E2 have a very small variance.

(ref:plot14env) Environmental variables with *low* variation and negative correlation. (a) E1 changing over time. (b) E2 changing over time.
```{r plot14env, fig.cap='(ref:plot14env)', fig.align="center", fig.height=6, fig.width=12}
library(MASS)
sample_size <- 15                                       
sample_meanvector <- c(mean(E1_series),  mean(E2_series))                                   
sample_covariance_matrix <- matrix(c(1, -1,
                                     -1, 1),
                                   ncol = 2)

# create bivariate normal distribution
set.seed(65354)
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
```



```{r}
# create communities
s <- 4
# community 1 - low variation in z1. 
z1_range <- 1.5
z2_range <- 50

# create reproducible community
set.seed(2456)
comm <- get_community(s, z1_range, z2_range)

# Get directional dreivatives
dir_1 <- get_directionals(comm, new_data, refs)


# community 2 - intermediate variation in z1. 
z1_range <- 25
z2_range <- 25
# b1_range <- 0.009
# b2_range <- 0.001

# create reproducible community
set.seed(2456)
comm <- get_community(s, z1_range, z2_range)

# Get directional dreivatives
dir_2 <- get_directionals(comm, new_data, refs)


# community 3 - large variation in z1. 
z1_range <- 50
z2_range <- 1.5
# b1_range <- 0.02
# b2_range <- 0.001

# create reproducible community
set.seed(2456)
comm <- get_community(s, z1_range, z2_range)

# Get directional dreivatives
dir_3 <- get_directionals(comm, new_data, refs)
```


### Response diversity calculation 

(ref:RDplot14) Directional derivatives and response diversity for the three different communities. 
a, b, c. Species directional derivatives over time. d, e, f. Response diversity measured as similarity-based diversity metric. g, h, i. Response diversity measured as divergence (sign sensitive). 
```{r RDplot14, fig.cap='(ref:RDplot14)', fig.align="center", fig.height=10, fig.width=16}
# community 1  - calculation RD

# from long to wide
rdiv_1 <- dir_1 %>%
  spread( sp, dir_deriv)
rdiv_1[is.na(rdiv_1)] <- 0

# actual calculation for only the same species used above
rdiv_1$rdiv<-apply(rdiv_1[,c(4:7)], 1, resp_div, sign_sens = F)
rdiv_1$sign<-apply(rdiv_1[,c(4:7)], 1, resp_div, sign_sens = T)
rdiv_1$Med<-mean(rdiv_1$rdiv)
rdiv_1$Med_sing<-mean(rdiv_1$sign)

# community 2  - calculation RD
# from long to wide
rdiv_2 <- dir_2 %>%
  spread( sp, dir_deriv)
rdiv_2[is.na(rdiv_2)] <- 0

# actual calculation for only the same species used above
rdiv_2$rdiv<-apply(rdiv_2[,c(4:7)], 1, resp_div, sign_sens = F)
rdiv_2$sign<-apply(rdiv_2[,c(4:7)], 1, resp_div, sign_sens = T)
rdiv_2$Med<-mean(rdiv_2$rdiv)
rdiv_2$Med_sing<-mean(rdiv_2$sign)
# community 3  - calculation RD
# from long to wide
rdiv_3 <- dir_3 %>%
  spread( sp, dir_deriv)
rdiv_3[is.na(rdiv_3)] <- 0

# actual calculation for only the same species used above
rdiv_3$rdiv<-apply(rdiv_3[,c(4:7)], 1, resp_div, sign_sens = F)
rdiv_3$sign<-apply(rdiv_3[,c(4:7)], 1, resp_div, sign_sens = T)
rdiv_3$Med<-mean(rdiv_3$rdiv)
rdiv_3$Med_sing<-mean(rdiv_3$sign)

plot_RD(dir_1, dir_2, dir_3, rdiv_1, rdiv_2, rdiv_3)

```

## High fluctuation in environmental change - positive correlation
Increasing variance in z1 and decreasing in z2 - positive correlation between env variables (E1 and E2) with *high* fluctuations

Same steps has before (3 communities with increasing diversity in z1 and z2), but there is positive correlation between E1 and E1


(ref:plot15env) Environmental variables with high variation and positive correlation. (a) E1 changing over time. (b) E2 changing over time.
```{r plot15env, fig.cap='(ref:plot15env)', fig.align="center", fig.height=6, fig.width=12}
library(MASS)
sample_size <- 15                                       
sample_meanvector <- c(mean(E1_series),  mean(E2_series))                                   
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
```


```{r}

# create communities
s <- 4
# community 1 - low variation in z1. 
z1_range <- 1.5
z2_range <- 50

# create reproducible community
set.seed(2456)
comm <- get_community(s, z1_range, z2_range)

# Get directional dreivatives
dir_1 <- get_directionals(comm, new_data, refs)

# community 2 - intermediate variation in z1. 
z1_range <- 25
z2_range <- 20
# b1_range <- 0.009
# b2_range <- 0.001

# create reproducible community
set.seed(2456)
comm <- get_community(s, z1_range, z2_range)

# Get directional dreivatives
dir_2 <- get_directionals(comm, new_data, refs)

# community 3 - large variation in z1. 
z1_range <- 50
z2_range <- 2
# b1_range <- 0.02
# b2_range <- 0.001

# create reproducible community
set.seed(2456)
comm <- get_community(s, z1_range, z2_range)

# Get directional dreivatives
dir_3 <- get_directionals(comm, new_data, refs)

```



### Response diversity calculation 

(ref:RDplot15) Directional derivatives and response diversity for the three different communities. 
a, b, c. Species directional derivatives over time. d, e, f. Response diversity measured as similarity-based diversity metric. g, h, i. Response diversity measured as divergence (sign sensitive). 
```{r RDplot15, fig.cap='(ref:RDplot15)', fig.align="center", fig.height=10, fig.width=16}
# community 1  - calculation RD

# from long to wide
rdiv_1 <- dir_1 %>%
  spread( sp, dir_deriv)
rdiv_1[is.na(rdiv_1)] <- 0

# actual calculation for only the same species used above
rdiv_1$rdiv<-apply(rdiv_1[,c(4:7)], 1, resp_div, sign_sens = F)
rdiv_1$sign<-apply(rdiv_1[,c(4:7)], 1, resp_div, sign_sens = T)
rdiv_1$Med<-mean(rdiv_1$rdiv)
rdiv_1$Med_sing<-mean(rdiv_1$sign)

# community 2  - calculation RD
# from long to wide
rdiv_2 <- dir_2 %>%
  spread( sp, dir_deriv)
rdiv_2[is.na(rdiv_2)] <- 0

# actual calculation for only the same species used above
rdiv_2$rdiv<-apply(rdiv_2[,c(4:7)], 1, resp_div, sign_sens = F)
rdiv_2$sign<-apply(rdiv_2[,c(4:7)], 1, resp_div, sign_sens = T)
rdiv_2$Med<-mean(rdiv_2$rdiv)
rdiv_2$Med_sing<-mean(rdiv_2$sign)

# community 3  - calculation RD
# from long to wide
rdiv_3 <- dir_3 %>%
  spread( sp, dir_deriv)
rdiv_3[is.na(rdiv_3)] <- 0

# actual calculation for only the same species used above
rdiv_3$rdiv<-apply(rdiv_3[,c(4:7)], 1, resp_div, sign_sens = F)
rdiv_3$sign<-apply(rdiv_3[,c(4:7)], 1, resp_div, sign_sens = T)
rdiv_3$Med<-mean(rdiv_3$rdiv)
rdiv_3$Med_sing<-mean(rdiv_3$sign)

plot_RD(dir_1, dir_2, dir_3, rdiv_1, rdiv_2, rdiv_3)
```


## Low fluctuation in environmental change - positive correlation 
Increasing variance in z1 and decreasing in 2 - positive correlation between env variables (E1 and E2) with *low* fluctuations.

We repeat the same steps as in the previous section, but this time E1 and E2 have a very small variance.

(ref:plot16env) Environmental variables with *low* variation and positive correlation. (a) E1 changing over time. (b) E2 changing over time.
```{r plot16env, fig.cap='(ref:plot16env)', fig.align="center", fig.height=6, fig.width=12}
library(MASS)
sample_size <- 15                                       
sample_meanvector <- c(mean(E1_series),  mean(E2_series))                                   
sample_covariance_matrix <- matrix(c(1, 1,
                                     1, 1),
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
```


```{r}

# create communities
s <- 4
# community 1 - low variation in z1. 
z1_range <- 1.5
z2_range <- 45

# create reproducible community
set.seed(2456)
comm <- get_community(s, z1_range, z2_range)

# Get directional dreivatives
dir_1 <- get_directionals(comm, new_data, refs)


# community 2 - intermediate variation in z1. 
z1_range <- 25
z2_range <- 20
# b1_range <- 0.009
# b2_range <- 0.001

# create reproducible community
set.seed(2456)
comm <- get_community(s, z1_range, z2_range)

# Get directional dreivatives
dir_2 <- get_directionals(comm, new_data, refs)


# community 3 - large variation in z1. 
z1_range <- 50
z2_range <- 2
# b1_range <- 0.02
# b2_range <- 0.001

# create reproducible community
set.seed(2456)
comm <- get_community(s, z1_range, z2_range)

# Get directional dreivatives
dir_3 <- get_directionals(comm, new_data, refs)

```



### Response diversity calculation 

(ref:RDplot16) Directional derivatives and response diversity for the three different communities. 
a, b, c. Species directional derivatives over time. d, e, f. Response diversity measured as similarity-based diversity metric. g, h, i. Response diversity measured as divergence (sign sensitive). 
```{r RDplot16, fig.cap='(ref:RDplot16)', fig.align="center", fig.height=10, fig.width=16}
# community 1  - calculation RD

# from long to wide
rdiv_1 <- dir_1 %>%
  spread( sp, dir_deriv)
rdiv_1[is.na(rdiv_1)] <- 0

# actual calculation for only the same species used above
rdiv_1$rdiv<-apply(rdiv_1[,c(4:7)], 1, resp_div, sign_sens = F)
rdiv_1$sign<-apply(rdiv_1[,c(4:7)], 1, resp_div, sign_sens = T)
rdiv_1$Med<-mean(rdiv_1$rdiv)
rdiv_1$Med_sing<-mean(rdiv_1$sign)

# community 2  - calculation RD
# from long to wide
rdiv_2 <- dir_2 %>%
  spread( sp, dir_deriv)
rdiv_2[is.na(rdiv_2)] <- 0

# actual calculation for only the same species used above
rdiv_2$rdiv<-apply(rdiv_2[,c(4:7)], 1, resp_div, sign_sens = F)
rdiv_2$sign<-apply(rdiv_2[,c(4:7)], 1, resp_div, sign_sens = T)
rdiv_2$Med<-mean(rdiv_2$rdiv)
rdiv_2$Med_sing<-mean(rdiv_2$sign)

# community 3  - calculation RD
# from long to wide
rdiv_3 <- dir_3 %>%
  spread( sp, dir_deriv)
rdiv_3[is.na(rdiv_3)] <- 0

# actual calculation for only the same species used above
rdiv_3$rdiv<-apply(rdiv_3[,c(4:7)], 1, resp_div, sign_sens = F)
rdiv_3$sign<-apply(rdiv_3[,c(4:7)], 1, resp_div, sign_sens = T)
rdiv_3$Med<-mean(rdiv_3$rdiv)
rdiv_3$Med_sing<-mean(rdiv_3$sign)

plot_RD(dir_1, dir_2, dir_3, rdiv_1, rdiv_2, rdiv_3)
```


## High fluctuation in environmental change - no correlation
Increasing variance in z1 and decreasing in z2 - no correlation between env variables (E1 and E2) with *high* fluctuations

Same steps has before (3 communities with increasing diversity in z1 and z2), but there is no correlation between E1 and E1


(ref:plot17env) Environmental variables with high variation and no correlation. (a) E1 changing over time. (b) E2 changing over time.
```{r plot17env, fig.cap='(ref:plot17env)', fig.align="center", fig.height=6, fig.width=12}
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
```


```{r}
# create communities
s <- 4
# community 1 - low variation in z1. 
z1_range <- 1.5
z2_range <- 45

# create reproducible community
set.seed(2456)
comm <- get_community(s, z1_range, z2_range)

# Get directional dreivatives
dir_1 <- get_directionals(comm, new_data, refs)

# community 2 - intermediate variation in z1. 
z1_range <- 25
z2_range <- 20
# b1_range <- 0.009
# b2_range <- 0.001

# create reproducible community
set.seed(2456)
comm <- get_community(s, z1_range, z2_range)

# Get directional derivatives
dir_2 <- get_directionals(comm, new_data, refs)


# community 3 - large variation in z1. 
z1_range <- 50
z2_range <- 1.5
# b1_range <- 0.02
# b2_range <- 0.001

# create reproducible community
set.seed(2456)
comm <- get_community(s, z1_range, z2_range)

# Get directional dreivatives
dir_3 <- get_directionals(comm, new_data, refs)

```

### Response diversity calculation 

(ref:RDplot17) Directional derivatives and response diversity for the three different communities. 
a, b, c. Species directional derivatives over time. d, e, f. Response diversity measured as similarity-based diversity metric. g, h, i. Response diversity measured as divergence (sign sensitive). 
```{r RDplot17, fig.cap='(ref:RDplot17)', fig.align="center", fig.height=10, fig.width=16}
# community 1  - calculation RD

# from long to wide
rdiv_1 <- dir_1 %>%
  spread( sp, dir_deriv)
rdiv_1[is.na(rdiv_1)] <- 0

# actual calculation for only the same species used above
rdiv_1$rdiv<-apply(rdiv_1[,c(4:7)], 1, resp_div, sign_sens = F)
rdiv_1$sign<-apply(rdiv_1[,c(4:7)], 1, resp_div, sign_sens = T)
rdiv_1$Med<-mean(rdiv_1$rdiv)
rdiv_1$Med_sing<-mean(rdiv_1$sign)

# community 2  - calculation RD
# from long to wide
rdiv_2 <- dir_2 %>%
  spread( sp, dir_deriv)
rdiv_2[is.na(rdiv_2)] <- 0

# actual calculation for only the same species used above
rdiv_2$rdiv<-apply(rdiv_2[,c(4:7)], 1, resp_div, sign_sens = F)
rdiv_2$sign<-apply(rdiv_2[,c(4:7)], 1, resp_div, sign_sens = T)
rdiv_2$Med<-mean(rdiv_2$rdiv)
rdiv_2$Med_sing<-mean(rdiv_2$sign)

# community 3  - calculation RD
# from long to wide
rdiv_3 <- dir_3 %>%
  spread( sp, dir_deriv)
rdiv_3[is.na(rdiv_3)] <- 0

# actual calculation for only the same species used above
rdiv_3$rdiv<-apply(rdiv_3[,c(4:7)], 1, resp_div, sign_sens = F)
rdiv_3$sign<-apply(rdiv_3[,c(4:7)], 1, resp_div, sign_sens = T)
rdiv_3$Med<-mean(rdiv_3$rdiv)
rdiv_3$Med_sing<-mean(rdiv_3$sign)

plot_RD(dir_1, dir_2, dir_3, rdiv_1, rdiv_2, rdiv_3)
```

## Low fluctuation in environmental change - no correlation
Increasing variance in z1 and decreasing in z2 - no correlation between env variables (E1 and E2) with *low* fluctuations

We repeat the same steps as in the previous section, but this time E1 and E2 have a very small variance.

(ref:plot18env) Environmental variables with *low* variation and negative correlation. (a) E1 changing over time. (b) E2 changing over time.
```{r plot18env, fig.cap='(ref:plot18env)', fig.align="center", fig.height=6, fig.width=12}
library(MASS)
sample_size <- 15                                       
sample_meanvector <- c(mean(E1_series),  mean(E2_series))                                   
sample_covariance_matrix <- matrix(c(1, 0.0,
                                     0.0, 1),
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
```


```{r}
# create communities
s <- 4
# community 1 - low variation in z1. 
z1_range <- 1.5
z2_range <- 45

# create reproducible community
set.seed(2456)
comm <- get_community(s, z1_range, z2_range)

# Get directional derivatives
dir_1 <- get_directionals(comm, new_data, refs)

# community 2 - intermediate variation in z1. 
z1_range <- 25
z2_range <- 20
# b1_range <- 0.009
# b2_range <- 0.001

# create reproducible community
set.seed(2456)
comm <- get_community(s, z1_range, z2_range)

# Get directional dreivatives
dir_2 <- get_directionals(comm, new_data, refs)

# community 3 - large variation in z1. 
z1_range <- 50
z2_range <- 2
# b1_range <- 0.02
# b2_range <- 0.001

# create reproducible community
set.seed(2456)
comm <- get_community(s, z1_range, z2_range)

# Get directional dreivatives
dir_3 <- get_directionals(comm, new_data, refs)
```




### Response diversity calculation 

(ref:RDplot18) Directional derivatives and response diversity for the three different communities. 
a, b, c. Species directional derivatives over time. d, e, f. Response diversity measured as similarity-based diversity metric. g, h, i. Response diversity measured as divergence (sign sensitive). 
```{r RDplot18, fig.cap='(ref:RDplot18)', fig.align="center", fig.height=10, fig.width=16}
# community 1  - calculation RD

# from long to wide
rdiv_1 <- dir_1 %>%
  spread( sp, dir_deriv)
rdiv_1[is.na(rdiv_1)] <- 0

# actual calculation for only the same species used above
rdiv_1$rdiv<-apply(rdiv_1[,c(4:7)], 1, resp_div, sign_sens = F)
rdiv_1$sign<-apply(rdiv_1[,c(4:7)], 1, resp_div, sign_sens = T)
rdiv_1$Med<-mean(rdiv_1$rdiv)
rdiv_1$Med_sing<-mean(rdiv_1$sign)

# community 2  - calculation RD
# from long to wide
rdiv_2 <- dir_2 %>%
  spread( sp, dir_deriv)
rdiv_2[is.na(rdiv_2)] <- 0

# actual calculation for only the same species used above
rdiv_2$rdiv<-apply(rdiv_2[,c(4:7)], 1, resp_div, sign_sens = F)
rdiv_2$sign<-apply(rdiv_2[,c(4:7)], 1, resp_div, sign_sens = T)
rdiv_2$Med<-mean(rdiv_2$rdiv)
rdiv_2$Med_sing<-mean(rdiv_2$sign)

# community 3  - calculation RD
# from long to wide
rdiv_3 <- dir_3 %>%
  spread( sp, dir_deriv)
rdiv_3[is.na(rdiv_3)] <- 0

# actual calculation for only the same species used above
rdiv_3$rdiv<-apply(rdiv_3[,c(4:7)], 1, resp_div, sign_sens = F)
rdiv_3$sign<-apply(rdiv_3[,c(4:7)], 1, resp_div, sign_sens = T)
rdiv_3$Med<-mean(rdiv_3$rdiv)
rdiv_3$Med_sing<-mean(rdiv_3$sign)

plot_RD(dir_1, dir_2, dir_3, rdiv_1, rdiv_2, rdiv_3)

```
