#Rarefaction script

#setting directory in your machine
setwd ("C:/Users/erick/Desktop/")

#loading required library

library (vegan)
#data input
asv_bact=read.csv("BroaMO_16S.csv",header=T,row.names="asv")
tax <- asv_bact[,13:19]
asv_bact <- asv_bact[,1:12]
t_asv_bact=t(asv_bact)

#rarefy
bac_rarefy<-rrarefy(t_asv_bact,min(rowSums(t_asv_bact)))
t_bac_rarefy=t(bac_rarefy)

dim(asv_bact);sum(asv_bact)
dim(bac_rarefy); sum(bac_rarefy)

#cut low abundant seqs
bact_cut<-t(bac_rarefy)
bact_cut<-bac_rarefy[,colSums(bac_rarefy) > 10] #test what is the best option here

dim(bac_rarefy); sum(bac_rarefy)
dim(bact_cut);sum(bact_cut)


bact_ra <- as.data.frame(bact_cut/rowSums(bact_cut))
rowSums(bact_ra)

dim(bact_ra);sum(bact_ra)

bact_ra=t(bact_ra)
bact_cut=t(bact_cut)


#Merge
x <- merge (bact_cut,tax,0,all.x = F)
z <- merge (bact_ra,tax,0,all.x = F)


write.csv (x,"BroaMO_16S_cut_rarefied+taxonomy.csv")
write.csv (z,"BroaMO_16S_ra_rarefied+taxonomy.csv")
