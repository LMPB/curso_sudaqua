#abundance variation per time
library(dplyr)
library(tidyr)
library(ggplot2)
setwd("C:/Users/erick/Desktop/")

asv_bact=read.csv("BroaMO_16S_ra_rarefied+taxonomy.csv",header=T,row.names="asv")
env=read.csv("BroaMO_env.csv")

groups = aggregate(asv_bact[1:12],by=list(category=asv_bact$Phylum), FUN = sum)
row.names(groups)=groups$category
groups=groups[,-1]

groups=t(groups)
groups=as.data.frame(groups)
groups[,25]<-env$date
colnames(groups)[25]<- "date"

nn=gather(as.data.frame(groups),key = "group", value = "value",-date) %>% 
  transform(group=factor(group,levels=c(
  "Abditibacteriota", "Acidobacteriota", "Actinobacteriota", "Armatimonadota", "Bacteroidota", "Bdellovibrionota",
"Chloroflexi", "Cyanobacteria", "Dependentiae", "Desulfobacterota", "FCPU426", "Firmicutes", "Fusobacteriota",
"Gemmatimonadota", "Hydrogenedentes", "Latescibacterota", "Margulisbacteria", "Myxococcota", "Patescibacteria",
"Planctomycetota", "Proteobacteria", "SAR324 clade(Marine group B)", "Spirochaetota", "Verrucomicrobiota"))) %>% 
  mutate(date = as.Date(date))


ggplot(nn,aes(x=date,y=value,fill=group))+
  geom_area()+
  scale_fill_manual(name="", values = c("green","blueviolet","firebrick","gainsboro","goldenrod",
                                        "mediumpurple","skyblue","cyan","sienna","violetred",
                                        "dodgerblue", "maroon4","seagreen","tomato","tan",
                                        "orange", "aliceblue","orchid","pink","red",
                                        "gray20", "gray50", "olivedrab","purple"))+ 
  theme_bw()+
  theme(axis.title.x=element_blank())+
  ylab(paste("Relative Abundances"))


