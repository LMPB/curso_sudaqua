#Betadispersion + nMDS
library(vegan)
library(dplyr)
library(tidyr)
library(cowplot)
library(gridExtra)
library(ggpubr)


setwd("C:/Users/erick/Desktop/")

asv_bact=read.csv("BroaMO_16S_ra_rarefied+taxonomy.csv",header=T,row.names="asv")
asv_bact=asv_bact[1:12]
asv_bact=t(asv_bact)

env=read.csv("BroaMO_env.csv", header=T, row.names= 1)
map = env[,11:14]

#Significance tests
dist_bray = vegdist(asv_bact, method="bray")
beta_dist_site=betadisper(dist_bray, map$trophic_state)
permutest(beta_dist_site, control=permControl (nperm=1000))
anova(beta_dist_site)
#TukeyHSD(beta_dist_site)
#plot(TukeyHSD(beta_dist_site), las = 1)

#extract axis' importance
valor = beta_dist_site$eig
valor2 = (beta_dist_site$eig/sum(beta_dist_site$eig))*100
head(valor2)

#AQUI
#extract axis' values
PCoA_score = scores(beta_dist_site$vectors[,1:2])
map_dist = merge(map, PCoA_score, by.x="row.names", by.y="row.names")

panel.a =
ggplot(map_dist2, aes(PCoA1, PCoA2, color = season)) + 
  geom_jitter(size=1.5) +  
  #scale_shape_manual(values=c(1,2,14,15,3,4,16,17,10)) + 
  scale_color_manual(name="Season", 
                     values=c("red2","green4", "navy","dodgerblue2","magenta2","chocolate3","purple"))+ 
  stat_ellipse(aes(group=season,color=season),size=0.75)+
  annotate("text", x = 0, y = -0.5, label = "PERMDISP = 0.080")+
  annotate("text", x = 0, y = -0.52, label = "PERMANOVA = 0.079")+
  theme_classic()+
  xlab("PCoA 1 (33.29%)")+
  ylab("PCoA 2 (17.18%)")

panel.b = 
ggplot(map_dist2, aes(PCoA1, PCoA2, color = trophic_state)) + 
  geom_jitter(size=1.5) +  
  #scale_shape_manual(values=c(1,2,14,15,3,4,16,17,10)) + 
  scale_color_manual(name="Trophic State", 
                     values=c("red2","green4", "navy","dodgerblue2","magenta2","chocolate3","purple"))+ 
  stat_ellipse(aes(group=trophic_state,color=trophic_state),size=0.75)+
  annotate("text", x = 0, y = -0.5, label = "PERMDISP = 0.093")+
  annotate("text", x = 0, y = -0.52, label = "PERMANOVA = 0.078")+
  theme_classic()+
  xlab("PCoA 1 (33.29%)")+
  ylab("PCoA 2 (17.18%)")

x = arrangeGrob(panel.a,panel.b,
                ncol = 1, nrow = 2,
                layout_matrix = rbind(c(1,2), c(1,2)))
as_ggplot(x)+
  draw_plot_label(label = c ("A","B"), size = 15,
                  x = c(0.4, 0.9), y = c(0.95,0.95))

#nMDS
ord_asv = metaMDS(asv_bact, dist="bray", k=2, trymax=200, autotransform=FALSE, noshare=FALSE, wascores=FALSE)
ord_asv$stress
nmds_score = scores(ord_asv)
map_dist = merge(map, nmds_score, by.x="row.names", by.y="row.names")

panel.a =
ggplot(map_dist, aes(NMDS1, NMDS2)) + 
  geom_jitter(aes(color = season), size=2) +  
  #scale_shape_manual(values=c(16,16,16)) + 
  scale_color_manual(name="Season", 
                     values=c("red2","green4", "navy","dodgerblue2","magenta2","chocolate3","purple"))+ 
  stat_ellipse(aes(group=season,color=season),size=0.75)+
  #geom_polygon(aes(group=season,fill=season),alpha=0.3)+
  geom_text(aes(label=Row.names),hjust=0, vjust=0, size=2)+
  annotate("text", x = 0, y = -0.5, label = "Stress = 0.052")+
  annotate("text", x = 0, y = -0.55, label = "PERMANOVA = 0.078")+
  theme_bw()+
  xlab("nMDS 1")+
  ylab("nMDS 2")

panel.b =
ggplot(map_dist, aes(NMDS1, NMDS2)) + 
  geom_jitter(aes(color = trophic_state), size=2) +  
  #scale_shape_manual(values=c(16,16,16)) + 
  scale_color_manual(name="Trophic State", 
                     values=c("red2","green4", "navy","dodgerblue2","magenta2","chocolate3","purple"))+ 
  stat_ellipse(aes(group=trophic_state,color=trophic_state),size=0.75)+
  #geom_polygon(aes(group=trophic_state,fill=trophic_state),alpha=0.3)+
  geom_text(aes(label=Row.names),hjust=0, vjust=0, size=2)+
  annotate("text", x = 0, y = -0.5, label = "Stress = 0.052")+
  annotate("text", x = 0, y = -0.55, label = "PERMANOVA = 0.079")+
  theme_bw()+
  xlab("nMDS 1")+
  ylab("nMDS 2")

x = arrangeGrob(panel.a,panel.b,
                ncol = 1, nrow = 2,
                layout_matrix = rbind(c(1,2), c(1,2)))
as_ggplot(x)+
  draw_plot_label(label = c ("A","B"), size = 15,
                  x = c(0.4, 0.9), y = c(0.95,0.95))
