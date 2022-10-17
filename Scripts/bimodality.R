library(dplyr)
library(tidyr)
library(ggplot2)
library(vegan)
library(cowplot)
library(gridExtra)
library(ggpubr)

setwd("C:/Users/erick/Desktop/")

asv_bact=read.csv("BroaMO_16S_ra_rarefied+taxonomy.csv",header=T,row.names="asv")
asv_bact=asv_bact[1:12]
asv_bact=t(asv_bact)
asv_bact_pa = decostand(asv_bact, "pa")
x=asv_bact_pa

a = print(ggplot(data=data.frame(apply(x, 1, sum)[which(apply(x, 2, sum) >0)]),aes(x=data.frame(apply(x, 2, sum)[which(apply(x, 2, sum) >0)])[,1]))+
            geom_histogram(fill="grey80",colour="black",binwidth=1)+
            ylim(c(0,ncol(x)))+
            ylab(paste("Frequency","(Number of ASVs)",sep="\n"))+
            xlab(paste("Persistence","(Number of times present)",sep="\n"))+
            ggtitle("")+
            coord_cartesian(ylim = c(0,220))+
            theme_bw())

table.out=as.data.frame(table(apply(x, 2, sum)[which(apply(x, 2, sum) >0)]))
colnames(table.out)= c("Persistency","Frequency")
print(table.out)
(ratio=(table.out[nrow(table.out),2])/(table.out[1,2]))

test.input=as.data.frame(table(apply(x, 2, sum)[which(apply(x, 2, sum) >0)]))
MOS.test=MOStest(as.numeric(table.out[,1]), 
                 as.numeric(table.out[,2]),
                 family=gaussian(link = "identity"), 
                 maxit = 200)
print(MOS.test)
capture.output(MOS.test,file=paste("bimodality_test.txt",sep="_")) #Saves table with stats

asv_bact=t(asv_bact)
asv_bact <- as.data.frame (asv_bact) %>% 
  mutate (mean = rowMeans(.)) %>% 
  mutate (count = rowSums(.[,-13] != 0))


b = print(ggplot(asv_bact,aes(x=count, y=mean))+
            geom_jitter(width = 0.8,alpha=0.8)+
            ylab(paste("Mean Relative Abundances"))+
            xlab(paste("Persistence","(Number of times present)",sep="\n"))+
            ggtitle("")+
            theme_bw())

x = arrangeGrob(a,b,
                ncol = 2, nrow = 1,
                layout_matrix = rbind(c(1,2)))
as_ggplot(x)+
  draw_plot_label(label = c ("A","B","*"), size = 15,
                  x = c(0, 0.52,0.45), y = c(1, 1,0.5))
