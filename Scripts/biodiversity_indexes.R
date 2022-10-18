#diversity indexes

#loading required libraries
library(BiodiversityR)
library(ggplot2)
library(dplyr)
library(tidyr)

#set your directory
setwd("C:/Users/erick/Desktop/")

#data input
asv_bact=read.csv("BroaMO_16S_ra_rarefied+taxonomy.csv",header=T,row.names="asv")
env= read.csv("Data/BroaMO_env.csv", header=T, row.names = 'id')

#Diversity Indexes
x=asv_bact[1:12]
x=t(x)
## "each site" for alfa-diversity; "pooled" for gamma-diversity
x1=diversityresult(x, y=NULL, level = NULL, factor = NULL,index = "richness", method = "each site", sortit = FALSE, digits = 5)
x2=diversityresult(x, y=NULL, level = NULL, factor = NULL,index = "Shannon", method = "each site", sortit = FALSE, digits = 5)
x3=diversityresult(x, y=NULL, level = NULL, factor = NULL,index = "Jevenness", method = "each site", sortit = FALSE, digits = 5)
x4=diversityresult(x, y=NULL, level = NULL, factor = NULL,index = "Simpson", method = "each site", sortit = FALSE, digits = 5)
x5=diversityresult(x, y=NULL, level = NULL, factor = NULL,index = "Eevenness", method = "each site", sortit = FALSE, digits = 5)

result=cbind(x1,x2,x3,x4,x5)
result[,1]=row.names(x1)
result[,2:6]=cbind(x1,x2,x3,x4,x5)
colnames(result)=c("site","richness","Shannon","Jevenness","Simpson","Eevenness")

write.table(result, "16s_diversity_results.txt")
result[,7]<-env$weather
colnames(result)[7]<- "weather"
result[,8]<-env$season
colnames(result)[8]<- "season"

#make de boxplot
nn=gather(as.data.frame(result),key = "index", value = "value",-site, -weather, -season) %>% 
  transform(index=factor(index,levels=c("richness","Shannon","Jevenness","Simpson","Eevenness")))

nn %>% 
  ggplot(aes(x=index,y=value))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width = 0.5, aes(color = season))+
  #facet_grid(~index, scales = "free")+
  facet_wrap(~index, ncol = 5, scales = "free")+
  theme_bw()+
  #theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(axis.title.x=element_blank())+
  theme(axis.text.x =element_blank())+
  theme(axis.ticks.x =element_blank())+
  theme(axis.title.y = element_blank())
