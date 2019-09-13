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
## specific richness south america freshwater
##
##########################################################################

## Libraries
library(countrycode)

## functions
source("00-scripts/step5/figures/functions.R")

## load data
dc = load_data("09-all_descripteurs/total_data_genetic_diversity_with_all_descripteurs.tsv")
dcT=dc[[1]];dcM=dc[[2]];dcF=dc[[3]];dcM.n1=dc[[4]];dcF.n1=dc[[5]];dcM.1=dc[[6]];dcF.1=dc[[7]]
richness=read.csv("01-infos/equalarea_id_coordsCA_FWRS_MR_RS.csv")

## FRESHWATER
Fr=merge(dcF.1,richness,by ="cell", all.x=F)
Fr=Fr[which(Fr$RS_FW_FISH !=0),]
Fr_wCoords=Fr[,c(2,3,4,6,7,27,10)] #with coordinates
coordsF=cbind(Fr_wCoords$x,Fr_wCoords$y)
#Fr$RS_FW_FISH=scale(Fr$RS_FW_FISH)
Fr=Fr[,c(1,4,6,7,27,10)]
Fr=na.omit(Fr)
names(Fr)=c("cell","GD_mean","nb_species","nb_indv_mean","richness","ISO3")

### add region continent information for FRESHWATER
north_america=c("USA","MEX","CAN")
regions=countrycode(Fr$ISO3,"iso3c","continent")
regions[which(Fr$ISO3 %in% north_america)]="North_america"
regions[which(regions=="Americas")]="South_america"
regions[which(is.na(regions))]="Antarctica"
Fr.r=cbind(Fr,regions)

summary(Fr.r[which(Fr.r$regions=="South_america"),]$richness)
