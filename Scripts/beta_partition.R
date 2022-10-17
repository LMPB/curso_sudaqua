#Beta-diversity Partitioning
library(betapart)
library(ggplot2)
library(vegan)
library(dplyr)
library(tidyr)
library(cowplot)
library(gridExtra)
library(ggpubr)

setwd("C:/Users/erick/Desktop/")


asv_bact=read.csv("BroaMO_16S_cut_rarefied+taxonomy.csv",header=T,row.names="asv")
asv_bact=asv_bact[1:12]
asv_bact=t(asv_bact)
asv_bact_pa = decostand(asv_bact, "pa")


#Beta-diversity partition
#The traditional betadiversity partition uses presence-absence data
#We also used the betapair abund using Bray-Curtis
beta.pair=beta.pair(asv_bact_pa, index.family = "jaccard")

nestedness <- data.frame(t(combn(rownames(asv_bact_pa),2)), as.numeric(beta.pair$beta.jne))
names(nestedness) <- c("S1", "S2", "nestedness")
turnover <- data.frame(t(combn(rownames(asv_bact_pa),2)), as.numeric(beta.pair$beta.jtu))
names(turnover) <- c("S1", "S2", "turnover")
dissimilarity <- data.frame(t(combn(rownames(asv_bact_pa),2)), as.numeric(beta.pair$beta.jac))
names(dissimilarity) <- c("S1", "S2", "dissimilarity")

#Boxplot
total <- merge (dissimilarity,turnover,by = c("S1","S2"), all = T)
complete <- merge (total,nestedness,by = c("S1","S2"), all = T)

panel.a = 
complete %>% 
  gather(key = "betadiv", value = "value", -S1, -S2) %>% 
  transform(betadiv=factor(betadiv,levels=c("dissimilarity","turnover","nestedness"))) %>% 
  ggplot(aes(x=betadiv,y=value, fill=betadiv))+
  geom_boxplot(outlier.shape = NA, width = 0.5,size=0.5, fill = "white")+
  geom_jitter(width = 0.2,alpha=0.2)+
  scale_fill_manual(name="", labels=c("dissimilarity","turnover","nestedness"), 
                    values = c("Grey90","Grey70","Grey50"))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(axis.title.x=element_blank())+
  theme(axis.title.y=element_blank())+
  theme(legend.position = "none")+
  ylim (0,1)


beta.pair.abund=beta.pair.abund(asv_bact, index.family = "bray")
abund.variation <- data.frame(t(combn(rownames(asv_bact),2)), as.numeric(beta.pair.abund$beta.bray.bal))
names(abund.variation) <- c("S1", "S2", "variation")
abund.gradient <- data.frame(t(combn(rownames(asv_bact),2)), as.numeric(beta.pair.abund$beta.bray.gra))
names(abund.gradient) <- c("S1", "S2", "gradient")
abund.dissimilarity <- data.frame(t(combn(rownames(asv_bact),2)), as.numeric(beta.pair.abund$beta.bray))
names(abund.dissimilarity) <- c("S1", "S2", "dissimilarity")


total_abund <- merge (abund.variation,abund.gradient,by = c("S1","S2"), all = T)
complete_abund <- merge (total_abund,abund.dissimilarity,by = c("S1","S2"), all = T)

panel.b =
complete_abund %>% 
  gather(key = "betadiv", value = "value", -S1, -S2) %>% 
  transform(betadiv=factor(betadiv,levels=c("dissimilarity","variation","gradient"))) %>% 
  ggplot(aes(x=betadiv,y=value, fill=betadiv))+
  geom_boxplot(outlier.shape = NA, width = 0.5,size=0.5, fill = "white")+
  geom_jitter(width = 0.2,alpha=0.2)+
  scale_fill_manual(name="", labels=c("dissimilarity","variation","gradient"), 
                    values = c("Grey90","Grey70","Grey50"))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(axis.title.x=element_blank())+
  theme(axis.title.y=element_blank())+
  theme(legend.position = "none")+
  ylim (0,1)


x = arrangeGrob(panel.a,panel.b,
                ncol = 1, nrow = 2,
                layout_matrix = rbind(c(1,2), c(1,2)))
as_ggplot(x)+
  draw_plot_label(label = c ("A","B"), size = 15,
                  x = c(0.4, 0.9), y = c(0.95,0.95))

