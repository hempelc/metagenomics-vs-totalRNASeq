library(gclus)
library(vegan)
library(dplyr)
library(data.table)
library(ellipse)
library(ape)
library(gdata)
library(ade4)
library(zCompositions)

#https://journals.plos.org/plosone/article/file?id=10.1371/journal.pone.0101238&type=printable
#https://help.xlstat.com/s/article/which-statistical-test-should-you-use?language=en_US

threshold<-0.0005 # discard entries with less than x% abundance #if veg is used need a threshold for veg? if not = 1
setwd("/home/gs/Schreibtisch/Soil_water")
############################################################
#####reading csv and transforming them into correct format##
############################################################
spe <- fread("soil_water_combined.csv_filtered0.000005info.csv",dec=",",data.table = F)
row.names(spe)<-spe$seqid
spe$V1<-NULL
spe_info<-spe[,1:16]
#samples by row for cmultrepl and clr
d.n0<-spe[,17:57] #soil
#d.n0<-spe[,58:98] #water
d.n0 = d.n0[ rowSums(d.n0)!=0, ]
d.n0<-t(d.n0)
d.n0 <- cmultRepl(d.n0, method="CZM", label=0)            # Chris - Replacing zeroes is useful, see package zCompositions
d.clr <- apply(d.n0, 1, function(x) log(x) - mean(log(x)))# perform clr on full dataset / centered log ratio
d.clrt<-t(d.clr)
spe<-as.data.frame(d.clrt)
row.names(spe)<-gsub(".x", "", rownames(spe))  # water was marked with samplename.y and soil with samplename.x previously
spe_soil_dist<-vegdist(spe,"euc") # euclidean distance on clr processed data = aitchison distance
