library(mosaic)
library(tidyverse)
View(KidsFeet)
KidsFeet %>% group_by(sex)
View(KidsFeet %>% group_by(sex))
  KidsFeet %>% group_by(sex) %>% 
  summarise(min=minlength),
  Q1=quantile(length,0.25),
  median=median(length),
  Q3=quantile(length,0.75),
  max=max(length))

KidsFeet %>% 
  group_by(sex) %>% 
  summarise(min=min(length))
  
