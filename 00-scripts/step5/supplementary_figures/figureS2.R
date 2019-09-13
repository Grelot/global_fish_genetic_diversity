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
## Supplementary Materials
## Figure S2 : Global distribution of higher and lower percentiles
## of genetic diversity: (a) and (c) 10th percentile
## and (b) and (d) 90th percentile
## of the distribution of genetic diversity
## for marine (a-b) and freshwater (c-d) fish species.
##
##########################################################################
## Libraries
lib_vect <-c("raster","sp","rgdal","rgeos","plyr","gridExtra","maptools")
sapply(lib_vect,library,character.only=TRUE)


##########################################################################
## functions
source("00-scripts/step5/figures/functions.R")

## from a grid and table of descriptors
## create a grid with Genetic diversity information for each cell of the grid
## attribute color "red" at mode="upper" or color "blue" at mode="upper"
grid_GD_mode <- function(gridW,dcW.1,dcW,mode) {
	gdW <- dcW[which(dcW$cell %in% dcW.1$cell),]$GD_mean
	gridW.1 <- gridW[which(gridW$IDcell %in% dcW.1$cell),]
	gridW.1@data$GD_mean <- gdW
	centroW <- gCentroid(gridW.1,byid=TRUE)
	if(mode=="upper") {             
		dfW=data.frame(lon=as.numeric(centroW@coords[,1]), lat=as.numeric(centroW@coords[,2]), col=rep("red",length(centroW@coords[,1])))
	} else {			
		dfW=data.frame(lon=as.numeric(centroW@coords[,1]), lat=as.numeric(centroW@coords[,2]), col=rep("blue",length(centroW@coords[,1])))
	}
	return(dfW)
}


## creates a grid with upper percentiles levels data for genetic diversity
grid_percentiles_upper <- function(dcV, dcV.1, gridV) {	
	upper90.1=quantile(dcV.1$GD_mean, 0.9)
	dcV.up=dcV.1[which(dcV.1$GD_mean >= upper90.1),]
	gridV.up=grid_GD_mode(gridV, dcV.up,dcV,"upper")
	return(gridV.up)
}


## creates a grid with lower percentiles levels data for genetic diversity
grid_percentiles_lower <- function(dcV, dcV.1, gridV) {
	lower10.1=quantile(dcV.1$GD_mean, 0.1)	
	dcV.low=dcV.1[which(dcV.1$GD_mean <= lower10.1),]
	gridV.low=grid_GD_mode(gridV, dcV.low,dcV,"lower")
	return(gridV.low)	
}


plot_GD_world <- function(gridW,worldcoast,biglakes,riverlakes,titleFig) {
    plot_world(worldcoast,biglakes,riverlakes)
    points(gridW$lon,gridW$lat,bg=as.character(gridW$col),pch=22,cex=0.65,col="black",lwd=0.2)
    box()
    mtext(titleFig,side=2,line=0.4,at=90,las=2)
}


##########################################################################
## load data
dc = load_data("09-all_descripteurs/total_data_genetic_diversity_with_all_descripteurs.tsv")
dcT=dc[[1]];dcM=dc[[2]];dcF=dc[[3]];dcM.n1=dc[[4]];dcF.n1=dc[[5]];dcM.1=dc[[6]];dcF.1=dc[[7]]
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
## only 90th percentile
gridM.up=grid_percentiles_upper(dcM,dcM.1,grid_Wgs84)
gridF.up=grid_percentiles_upper(dcF,dcF.1,grid_Wgs84)

## only 10th percentile
gridM.low=grid_percentiles_lower(dcM,dcM.1,grid_Wgs84)
gridF.low=grid_percentiles_lower(dcF,dcF.1,grid_Wgs84)


##########################################################################
## write pdf files
tiff(filename="10-figures/figureS2.tiff",width=22,height=16,units="cm",res=640,compression="lzw")

layout(mat=rbind(c(1,2),c(1,2),c(3,4),c(3,4)))
par(mar=c(2,2,2,2))
plot_GD_world(gridM.up,worldcoast,biglakes,riverlakes,"(a)")
plot_GD_world(gridM.low,worldcoast,biglakes,riverlakes,"(b)")
plot_GD_world(gridF.up,worldcoast,biglakes,riverlakes,"(c)")
plot_GD_world(gridF.low,worldcoast,biglakes,riverlakes,"(d)")

dev.off()
