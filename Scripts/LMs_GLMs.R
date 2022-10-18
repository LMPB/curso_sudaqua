##Generalized Linear Models

###### 1.1 - loading required libraries #####
{
library(ggplot2)
library(gtable)
library(scales)
library(gridExtra)
library(dplyr)
library(tidyr)
library(ggrepel)
library(grid)
library(ggpubr)
library(pastecs)
library(gtools)
}
###### 1.2 - External Functions #####
##Dsquare function: used to extract the R2 value from GLMs
#the external script is also available in the same github repository
source("00_Dsquare_GLM.R")

###### 2.1 - Data input ######
data <- read.table("16s_diversity_env_final.txt",header = TRUE)

###### 2.2 - Transform to long format
data_2 <- data %>%
  dplyr::select(-Row.names) %>%
  gather(key = variable, value = value,-site,-weather,-season,-trophic_state,-ph_scale,-date) %>%
  na.omit()

data_2$value<-as.numeric(data_2$value) # transform all values as numeric

#### 2.3 - Observing data distribution ####
data_2 %>% ggplot(aes(value))+
  facet_wrap(.~variable,scale="free")+
  geom_density(aes(fill="blue"))+
  theme_bw()+theme(legend.position='none')

#### 3.1.1 - Running Linear Models
m1 <- lm(Shannon ~ temp_water, data=data) 
summary(m1) #R2= 0.3071; p=0.03581
m2 <- lm(Shannon ~ ph, data=data) 
summary(m2)#N.S
m3 <- lm(Shannon ~ doc, data=data) 
summary(m3)#N.S
m4 <- lm(Shannon ~ tn, data=data) 
summary(m4)#N.S
m5 <- lm(Shannon ~ chla, data=data) 
summary(m5) #N.S.
m6 <- lm(Shannon ~ tsi, data=data) 
summary(m6) #R2=0.2096; p=0.1915

##### GLMs
m1.2<-glm(Shannon ~ temp_water, data=data, family=Gamma(link="log"))
summary(m1.2);Dsquared(m1.2) #0.36
m2.2<-glm(Shannon ~ tsi, data=data, family=Gamma(link="log"))
summary(m2.2);Dsquared(m2.2) #R2=0.2763124; p=0.08

#Multiple regressions
mult_model<-  lm (Shannon ~ temp_water*season, data=data) 
summary(mult_model)#N.S
mult_model<-  lm (Shannon ~ temp_water*weather, data=data) 
summary(mult_model)#N.S

###DONE###
