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
## functions used to generate figures
## Supplementary Figure 5. Sampling effect.
## (a,b) number of sequences per cell,
## (c,d) number of fish species per cell,
## (e,f) mean number of sequences per fish species per cell for marine species
## (a, c, e) and freshwater species (b, d,f) respectively. 
##
##
##########################################################################
## Libraries
lib_vect <-c("raster","plotrix","rgeos","rgdal","sp","maptools","shape","parallel","png","plyr")
sapply(lib_vect,library,character.only=TRUE)


##########################################################################
## functions
source("00-scripts/step5/figures/functions.R")

## convert number of sequences value into color value
color_N <- function(gd){
	gdc=gd
	gdc[which(gd<5)] = "< 5"	
	gdc[which(gd>=5 & gd<10)] = "5 - 10"
	gdc[which(gd>=10 & gd<20)] ="10 - 20"
	gdc[which(gd>=20 & gd<30)] ="20 - 30"	
	gdc[which(gd>=30)] =">30"
	return(gdc)
}

color_nspec <- function(gd){
	gdc=gd
	gdc[which(gd<3)] = "2 - 3"	
	gdc[which(gd>=3 & gd<5)] = "3 - 5"
	gdc[which(gd>=5 & gd<10)] ="5 - 10"	
	gdc[which(gd>=10)] =">10"
	return(gdc)
}

color_ind <- function(gd){
	gdc=gd
	gdc[which(gd<4)] = "2 - 4"	
	gdc[which(gd>=4 & gd<6)] = "4 - 6"
	gdc[which(gd>=6 & gd<8)] ="6 - 8"
	gdc[which(gd>=8 & gd<10)] ="8 - 10"	
	gdc[which(gd>=10)] =">10"
	return(gdc)
}


## from a grid and table of descriptors
## create a grid with number of sequences information for each cell of the grid
grid_numberofseq <- function(gridW,dcW.1,dcW) {
	indW=dcW[which(dcW$cell %in% dcW.1$cell),]$nb_indv_mean
	specW=dcW[which(dcW$cell %in% dcW.1$cell),]$nb_species
	nW=indW*specW
	gridW.1=gridW[which(gridW$IDcell %in% dcW.1$cell),]	
	centroW <- gCentroid(gridW.1,byid=TRUE)
	dfW=data.frame(lon=as.numeric(centroW@coords[,1]),
	 lat=as.numeric(centroW@coords[,2]),
	 numberOfSequences=color_N(nW),
	 numberOfSequences_color=color_N(nW),
	 numberOfSpecies=color_nspec(specW),
	 numberOfSpecies_color=color_nspec(specW),
	 meanNumberOfIndividualsBySpecies=color_ind(indW),
	 meanNumberOfIndividualsBySpecies_color=color_ind(indW),
	 nseq=nW,
	 nspec=specW,
	 nind=indW)
	dfW.ordered=dfW[order(dfW$nseq),]
	dfW.ordered$numberOfSequences=factor(dfW.ordered$numberOfSequences,
		levels=c("< 5","5 - 10", "10 - 20", "20 - 30",">30" ))
	dfW.ordered=dfW.ordered[order(dfW.ordered$nspec),]
	dfW.ordered$numberOfSpecies=factor(dfW.ordered$numberOfSpecies,
		levels=c("2 - 3","3 - 5", "5 - 10", ">10" ))
	dfW.ordered=dfW.ordered[order(dfW.ordered$nind),]
	dfW.ordered$meanNumberOfIndividualsBySpecies=factor(dfW.ordered$meanNumberOfIndividualsBySpecies,
		levels=c("2 - 4","4 - 6", "6 - 8", "8 - 10", ">10" ))
	dfW.ordered$numberOfSequences_color=dfW.ordered$numberOfSequences
	dfW.ordered$numberOfSpecies_color=dfW.ordered$numberOfSpecies
	dfW.ordered$meanNumberOfIndividualsBySpecies_color=dfW.ordered$meanNumberOfIndividualsBySpecies
	dfW.ordered=dfW.ordered[order(dfW.ordered$nseq),]
	levels(dfW.ordered$numberOfSequences_color) <- c("#FFFFFF","#FFF0F0","#FFB0B0","#FF7575","#FF0000")
	dfW.ordered=dfW.ordered[order(dfW.ordered$nspec),]
	levels(dfW.ordered$numberOfSpecies_color) <- c("#FFFFFF","#FFF0F0","#FF6F6F","#FF0000")
    dfW.ordered=dfW.ordered[order(dfW.ordered$nind),]
	levels(dfW.ordered$meanNumberOfIndividualsBySpecies_color) <- c("#FFFFFF","#FFF0F0","#FFB0B0","#FF7575","#FF0000")
	return(dfW.ordered)
}

## from a grid with number of sequences information for each cell
## and from a worldcoast shape polygon object
## and from a river and lakes shape polygon object
## create a plot map of number of sequences
plot_nseq_world <- function(gridW,colorColon,worldcoast,biglakes,riverlakes,titleFig) {
    plot_world(worldcoast,biglakes,riverlakes)
    points(gridW$lon,gridW$lat,bg=as.character(gridW[,colorColon]),pch=22,cex=0.65,col="black",lwd=0.2)
    box()
    mtext(titleFig,side=2,line=0.4,at=90,las=2)
}


##########################################################################
## load data
dc <- load_data("09-all_descripteurs/total_data_genetic_diversity_with_all_descripteurs.tsv")
dcT <- dc[[1]];dcM=dc[[2]];dcF=dc[[3]];dcM.n1=dc[[4]];dcF.n1=dc[[5]];dcM.1=dc[[6]];dcF.1=dc[[7]]
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
## filter wrong freshwater assignation
wrongFreshwaterCells=c(3304,4741,8453,9146,11400,11474,11548,12434,12730,11860)
dcF.1=dcF.1[which(!dcF.1$cell %in% wrongFreshwaterCells),]


##########################################################################
## load worldcoast shapefile
grid <- readOGR("01-infos/grid_equalarea200km","gridFish.b_260418")
worldcoast <- readOGR("01-infos/ne_50m_land",layer="ne_50m_land")

## load rivers and lake centerlines shapefile
riverlakes  <- readOGR("01-infos/ne_50m_rivers_lake_centerlines_scale_rank",layer="ne_50m_rivers_lake_centerlines_scale_rank")
grid_Wgs84 <- spTransform(grid, proj4string(worldcoast))
biglakes <- readShapePoly("01-infos/big_lakes/GSHHS_h_L2.shp",proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"))


##########################################################################
## add number of sequences to the grid cells
gridM.1.df=grid_numberofseq(grid_Wgs84,dcM.1,dcM)
gridF.1.df=grid_numberofseq(grid_Wgs84,dcF.1,dcF)


##########################################################################
## write pdf files
tiff(filename="10-figures/figureS5.tiff",width=27,height=24,units="cm",res=640,compression="lzw")
layout(mat=rbind(c(1,2),c(3,4),c(5,6)))
par(mar=c(1.5,2,2,1))

## number of sequences per cell
plot_nseq_world(gridM.1.df,4,worldcoast,biglakes,riverlakes,"(a)")
plot_nseq_world(gridF.1.df,4,worldcoast,biglakes,riverlakes,"(b)")
legend("bottomleft",legend=levels(gridM.1.df$numberOfSequences),col=levels(gridM.1.df$numberOfSequences_color),
       pch=15,cex=0.9,title="Number of sequences",pt.cex=1)
## number of species per cell
plot_nseq_world(gridM.1.df,6,worldcoast,biglakes,riverlakes,"(c)")
plot_nseq_world(gridF.1.df,6,worldcoast,biglakes,riverlakes,"(d)")
legend("bottomleft",legend=levels(gridM.1.df$numberOfSpecies),col=levels(gridM.1.df$numberOfSpecies_color),
       pch=15,cex=0.9,title="Number of species",pt.cex=1)
## number of individuals by species per cell
plot_nseq_world(gridM.1.df,8,worldcoast,biglakes,riverlakes,"(e)")
plot_nseq_world(gridF.1.df,8,worldcoast,biglakes,riverlakes,"(f)")
legg=legend("bottomleft",legend=levels(gridM.1.df$meanNumberOfIndividualsBySpecies),col=levels(gridM.1.df$meanNumberOfIndividualsBySpecies_color), pch=15,cex=0.9,title="Mean number of\nsequences by species",pt.cex=1,bty='n')
lxleft <- legg$rect[["left"]]
lytop <- legg$rect[["top"]]
lybottom <- lytop - legg$rect[["h"]]
lxright <- lxleft + legg$rect[["w"]]
rect(lxleft, lybottom-.1, lxright, lytop+6,col="white")
legend("bottomleft",legend=levels(gridM.1.df$meanNumberOfIndividualsBySpecies),col=levels(gridM.1.df$meanNumberOfIndividualsBySpecies_color), pch=15,cex=0.9,title="Mean number of\nsequences by species",pt.cex=1,bty='n')

dev.off()
