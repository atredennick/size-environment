library(tidyverse)
library(ggthemes)
library(ggalt)
library(mvtnorm)


### Negatively covarying vital rates
ys <- rmvnorm(n = 15, mean = rep(0, 2), matrix(data = c(1,-0.9,-0.9,1),
                                               ncol = 2, nrow = 2))
df <- as.data.frame(ys) %>%
  mutate(Year = as.integer(1:nrow(ys))) %>%
  gather(group, value, -Year)

ggplot(df, aes(Year,value,color=group))+
  geom_xspline(spline_shape = 1, size = 2.5)+
  scale_color_brewer(type = "qual")+
  xlab(NULL) + ylab(NULL)+
  guides(color = FALSE)+
  theme_hc()+
  theme(axis.ticks = element_blank(),
        axis.text = element_blank())



### Reduced variance of large plants
ys <- rmvnorm(n = 15, mean = rep(0, 2), matrix(data = c(1,0,0,1*0.01),
                                               ncol = 2, nrow = 2))
df <- as.data.frame(ys) %>%
  mutate(Year = as.integer(1:nrow(ys))) %>%
  gather(group, value, -Year)

ggplot(df, aes(Year,value,color=group))+
  geom_xspline(spline_shape = 1, size = 2.5)+
  scale_color_brewer(type = "qual")+
  xlab(NULL) + ylab(NULL)+
  guides(color = FALSE)+
  theme_hc()+
  theme(axis.ticks = element_blank(),
        axis.text = element_blank())

