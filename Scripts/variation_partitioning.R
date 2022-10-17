#Beta-diversity Partitioning

library(dplyr)
library(tidyr)
library(ggplot2)
library(betapart)
library(vegan)


setwd("C:/Users/erick/Desktop/")

asv_bact=read.csv("BroaMO_16S_ra_rarefied+taxonomy.csv",header=T,row.names="asv")
asv_bact=asv_bact[1:12]
asv_bact=t(asv_bact)
beta.pair.abund=beta.pair.abund(asv_bact, index.family = "bray")

env=read.csv("BroaMO_env.csv", header=T, row.names= 1)
denv=decostand(env[2:10], na.rm=T, method="standardize")
denv[,5]<-env[,6]

#Forward selection for each dataset
#Env variables selection
rda2env <- capscale(beta.pair.abund$beta.bray ~ ., data=denv,distance='jaccard')
rda1env <- capscale(beta.pair.abund$beta.bray ~ 1, data=denv, distance='jaccard')
step.forward <- ordiR2step(rda1env,scope=formula(rda2env), direction="forward", perm.max=200,pstep=999)
step.forward$anova
names(denv)

temp = denv[,c(2)]
secchi = denv[,c(3)]

var.part <- varpart(beta.pair.abund$beta.bray,temp,secchi)
plot(var.part, digits=2, Xnames = c('Temperature', 'Secchi'))
var.part
showvarparts(2)

#Test of all testable fractions
#test of fractions [a+b+c]
all.pars <-cbind(temp,secchi)
anova.cca(dbrda(beta.pair.abund$beta.bray~as.matrix(all.pars)),step=1000)
#test of fractions [a+b] => Environmental fraction
anova.cca(dbrda(beta.pair.abund$beta.bray~as.matrix(temp)),step=1000)
#test of fractions [b+c] => Spatial fraction
anova.cca(dbrda(beta.pair.abund$beta.bray~as.matrix(secchi)),step=1000)
#test of fraction [a]
anova.cca(dbrda(beta.pair.abund$beta.bray~as.matrix(temp)+ Condition(as.matrix(secchi)),sqrt.dist = T),step=1000)
#test of fraction [c]
anova.cca(dbrda(beta.pair.abund$beta.bray~as.matrix(secchi)+ Condition(as.matrix(temp)),sqrt.dist = T),step=1000)
