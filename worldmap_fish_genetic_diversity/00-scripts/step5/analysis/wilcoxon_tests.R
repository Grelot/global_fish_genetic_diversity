##########################################################################
## Codes for the paper:
## Global determinants of freshwater and marine fish genetic diversity
## Authors :
## Stephanie Manel, Pierre-Edouard Guerin, David Mouillot,
## Simon Blanchet, Laure Velez, Camille Albouy, Loic Pellissier
## 
## Montpellier, 2017-2019
## Submited to Nature communications, 2019
##
## wilcoxon test
## determination mediane 
##
##########################################################################
## Libraries
library(plyr)

##########################################################################
## load data
## read the table
dcT=read.table("09-all_descripteurs/total_data_genetic_diversity_with_all_descripteurs.tsv",header=T)
## keep entries with GD!=0
dcT=dcT[which(dcT$GD_mean != 0),]
## keep entries with more than 2 species
#cluster marine
dcM<-dcT[which(dcT$is_marine == 1),]
##cluster freshwater
dcF<-dcT[which(dcT$is_marine == 0),]
#remove cells with only one species
dcM.1<-dcM[which(dcM$nb_species != 1),]
dcF.1<-dcF[which(dcF$nb_species != 1),]
## remove marine data with no SST information
dcM.1=dcM.1[which(!is.na(dcM.1$SST)),]
## remove marine data with no CloMeanVal information
dcM.1=dcM.1[which(!is.na(dcM.1$cloMeanVal)),]
## remove marine data with no bathymetry information
dcM.1=dcM.1[which(!is.na(dcM.1$bathyVal)),]
## remove freshwater data with no mean regional temperature information
dcF.1=dcF.1[which(!is.na(dcF.1$SST)),]
## remove freshwater data with no bathyval information
dcF.1=dcF.1[which(!is.na(dcF.1$bathyVal)),]

## load richness data
richness=read.csv("01-infos/equalarea_id_coordsCA_FWRS_MR_RS.csv")
Fr=merge(dcF.1,richness,by ="cell", all.x=F)
Fr=Fr[which(Fr$RS_FW_FISH !=0),]
Fr_wCoords=Fr[,c(2,3,4,6,7,27)] #with coordinates
Fr=Fr[,c(1,4,6,7,27)]
Fr=na.omit(Fr)
names(Fr)=c("cell","GD_mean","nb_species","nb_indv_mean","richness")
Mr=merge(dcM.1,richness,by ="cell", all.x=F)
Mr=Mr[which(Mr$RS_MR_FISH !=0),]
Mr_wCoords=Mr[,c(2,3,4,6,7,28)] #with coordinates
Mr=Mr[,c(1,4,6,7,28)]
Mr=na.omit(Mr)
names(Mr)=c("cell","GD_mean","nb_species","nb_indv_mean","richness")
Tr=rbind(Mr,Fr)




##########################################################################
## Estimation of the median
summary(dcM.1)
summary(dcF.1)
summary(Mr$richness)
summary(Fr$richness)

boxplot(dcM.1$GD_mean,dcF.1$GD_mean)
hist(dcM.1$GD_mean)
## test de Wilcoxon Genetic diversity
wilcox.test(dcM.1$GD_mean,dcF.1$GD_mean) # W =69070, p-value = 1.101e-09
## Wilcoxon Specific richness # W = 79143, p-value = 0.7306
wilcox.test(Mr$richness, Fr$richness)