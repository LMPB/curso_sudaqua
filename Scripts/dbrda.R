#dbRDA
library(vegan)

setwd("C:/Users/erick/Desktop/")

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

dbRDA=capscale(asv_bact ~ temp_air+temp_water+secchi+ph+doc+tn+chla, denv, dist="bray")
#extract axis' importance and significance
anova(dbRDA)
valor = (dbRDA$CCA$eig/sum(dbRDA$CCA$eig))*100
head(valor)


#plot dbRDA
plot(dbRDA, type="n", xlim = c(-2, 2), ylim = c(-3, 3), xlab="DIM1 (36.48%)", ylab="DIM2 (23.81%)")
points(dbRDA, pch=16, col="black", bg="black", cex=0.8)
text(dbRDA, dis="bp", col="red3", cex=1, font = 3) #if "error in seq_len(nrow(pts)) argument must be coercible to non-negative integer", change 'cn' by 'bp'
text(dbRDA,"sites", col="black", cex=0.6, pos=1)
text("p=0.087", x=-3, y=-2)
