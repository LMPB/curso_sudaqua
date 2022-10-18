#Mantel test script

#loading required libraries
library(vegan)
library(dplyr)
library(tidyr)
library(ggplot2)

#set directory in your machine
setwd("C:/Users/erick/Desktop/")

#data input
asv_bact=read.csv("BroaMO_16S_ra_rarefied+taxonomy.csv",header=T,row.names="asv")
asv_bact=asv_bact[1:12]
asv_bact=t(asv_bact)

env=read.csv("BroaMO_env.csv", header=T, row.names= 1)
denv=decostand(env[2:10], na.rm=T, method="standardize")
denv[,5]<-env[,6]

rdasenv <- rda(asv_bact ~ temp_air+temp_water+secchi+ph+doc+tn+chla,data=denv)
vif.cca(rdasenv)#colinearidade: recomendável tirar valores maiores que 10 das análises

#Parsimonious subsets of explanatory variables (based on #forward selection)
names(denv)
denv <- denv[, c(1,2,3,5,6,7,8)]#seleciona as amostras que passaram pela análise de parsimônia e colinearidade
names(denv)


#Biologic Distance matrices
bac.dist <- vegdist(asv_bact, "bray") #para diversidade melhor usar Bray

#Environmental Distance matrices
#Mantel test: Abundance vs Environment 
env.dist <- vegdist(denv[,7], "euclidean") #aqui eu escolho uma variável ambiental da minha tabela para fazer a matriz de dissimilaridade
mantel(bac.dist, env.dist, method="spear") #mantel abundancias vs altitude


#Plot graph
mantel=read.table(file = "mantel.txt", header=TRUE, sep="\t", dec=".", na.strings="NA") #input mantel results for plotting

mantel %>% 
  gather(key = "type", value = "valor", -variable) %>% 
  transform(variable=factor(variable,levels=c("temp_water","secchi","chla"))) %>% 
  ggplot(aes(x=type,y=valor, group=variable, color=variable))+
  geom_point(size = 2)+
  geom_line()+
  scale_color_manual(name="", 
                     values = c("red2","cyan1","purple4","grey50", "green4","gold","purple","green"))+
  theme_bw()+
  scale_x_discrete(labels=c("Time"))+
  #scale_x_discrete(labels=c("Bimodal","Gamma","Other"))+
  xlab("")+
  ylab("")

##Heatmap
mantel %>% 
  gather(key = "type", value = "valor", -variable) %>% 
  #transform(type=factor(type,levels=c("Core","Non.core","Satellite"))) %>% 
  #transform(type=factor(type,levels=c("bimodal","gamma","other"))) %>% 
  transform(variable=factor(variable,levels=c("temp_water","secchi","chla"))) %>% 
  ggplot(aes(x=type,y=variable)) +
  geom_tile(aes(fill = valor))+
  theme_bw()+
  #scale_x_discrete(labels=c("Core","Non-core","Satellite"))+
  scale_x_discrete(labels=c("16S","18S"))+
  scale_fill_gradient(low = "white", high = "Navy")+
  xlab("")+
  ylab("")

