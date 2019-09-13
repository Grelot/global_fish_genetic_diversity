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
## Figure 1: Map of the global distribution of genetic diversity for marine species
## Average number of mutations per base pair for CO1 across
## (a) saltwater species
## (b) freshwater species. 
## The color gradient represent the relative genetic diversity :
## the reddest square area have the highest genetic diversity.
## The bluest square area have the lowest genetic diversity.
##
##
##########################################################################
## Libraries
lib_vect <-c("raster","plotrix","rgeos","rgdal","sp","maptools","shape","parallel","png","plyr")
sapply(lib_vect,library,character.only=TRUE)

##########################################################################
## functions
source("00-scripts/step5/figures/functions.R")

## convert genetic diversity value into color value
color_GD <- function(gd){
	gdc=gd
	gdc[which(gd<0.001)] = "< 0.001"	
	gdc[which(gd>=0.001 & gd<0.002)] = "0.001 - 0.002"
	gdc[which(gd>=0.002 & gd<0.003)] =" 0.002 - 0.003"
	gdc[which(gd>=0.003 & gd<0.005)] ="0.003 - 0.005"
	gdc[which(gd>=0.005 & gd<0.01)] ="0.005 - 0.01"
	gdc[which(gd>=0.01 & gd<0.02)] ="0.01 - 0.025"
	gdc[which(gd>0.02)] =">0.025"
	return(gdc)
}

## convert count cells value into color value
color_countcells <- function(cc){
	cdc=cc
	cdc[which(cc<10)] = "< 10"	
	cdc[which(cc>=10 & cc<20)] = "10 - 20"
	cdc[which(cc>=20 & cc<40)] ="20 - 40"
	cdc[which(cc>=40 & cc<60)] ="40 - 60"
	cdc[which(cc>60)] =">60"
	return(cdc)

}

## from a grid and table of descriptors
## create a grid with Genetic diversity information for each cell of the grid
grid_GD <- function(gridW,dcW.1,dcW) {
	gdW <- dcW[which(dcW$cell %in% dcW.1$cell),]$GD_mean
	gridW.1 <- gridW[which(gridW$IDcell %in% dcW.1$cell),]
	gridW.1@data$GD_mean <- gdW
	centroW <- gCentroid(gridW.1,byid=TRUE)
	dfW=data.frame(lon=as.numeric(centroW@coords[,1]), lat=as.numeric(centroW@coords[,2]), GD=color_GD(gdW),gdmean=gdW,col=color_GD(gdW))
	dfW.order=dfW[order(dfW$gdmean),]
	dfW.order$GD=factor(dfW.order$GD, levels=c("< 0.001","0.001 - 0.002", " 0.002 - 0.003", "0.003 - 0.005","0.005 - 0.01","0.01 - 0.025",">0.025" ))
	dfW.order$col=dfW.order$GD
	levels(dfW.order$col) <- c("#3333A2","#3333FF","#33CBFF","#33FFFF","#FFDF33","#FFA333","#FF3333")
	return(dfW.order)
}

## from a grid with genetic diversity information for each cell
## and from a worldcoast shape polygon object
## and from a river and lakes shape polygon object
## create a plot map of genetic diversity
plot_GD_world <- function(gridW,worldcoast,biglakes,riverlakes,titleFig) {
    plot_world(worldcoast,biglakes,riverlakes)
    points(gridW$lon,gridW$lat,bg=as.character(gridW$col),pch=22,cex=0.64,col="black",lwd=0.2)
    box()
    mtext(titleFig,side=2,line=0.4,at=96,las=2,cex=1.2)
}

## convert coordinates into epsg:4326 and aggregate by latband(10degree)
## and gives the mean of the genetic diversity per species per cell per latband
latbandGD <- function(dcW.1, dcW,gridW) {	
	gridW.1=gridW[which(gridW$IDcell %in% dcW.1$cell),]
	dcW=dcW[which(dcW$cell %in% dcW.1$cell),]
	centroW <- gCentroid(gridW.1,byid=TRUE)	
	spW=SpatialPoints(data.frame(as.numeric(centroW@coords[,1]),as.numeric(centroW@coords[,2])))
	proj4string(spW)=proj4string(gridW)
	spW.g=spTransform(spW, CRS("+init=epsg:4326"))
	lat_W=matrix(spW.g@coords[,2])
	latband_W=apply(lat_W, 1, FUN=function(x) coords_to_latband(x))
	dcW.latband=data.frame(latband=latband_W,GD=dcW$GD_mean)
	countCellsbyLatband=as.vector(table(dcW.latband$latband))
	latband.cc=as.integer(names(table(dcW.latband$latband)))
	dcW.lb.cc=data.frame(latband=latband.cc, countCells=color_countcells(countCellsbyLatband), col=color_countcells(countCellsbyLatband))
	dcW.lb.mn=aggregate(GD ~ latband, FUN = mean, dcW.latband)
	names(dcW.lb.mn)=c("latband","GD_band_mean")
	dcW.lb.ci=aggregate(GD ~ latband, FUN = function(x) sd(x,na.rm=TRUE)/sqrt(length(x)), dcW.latband)
    names(dcW.lb.ci)=c("latband","GD_band_ci")
	dcW.mnci=merge(dcW.lb.mn,dcW.lb.ci,by="latband",all=T)	
  	latband=data.frame(latband=seq(from=-90,to=90,by=10))
  	dcW.alb=merge(latband, dcW.mnci , by="latband",all=T)
	dcW.albc=merge(dcW.alb,dcW.lb.cc, by="latband",all=T)
	if(any(is.na(dcW.albc$GD_band_mean))){dcW.albc[which(is.na(dcW.albc$GD_band_mean)),]$GD_band_mean = 0 }
	if(any(is.na(dcW.albc$GD_band_ci))){dcW.albc[which(is.na(dcW.albc$GD_band_ci)),]$GD_band_ci = 0 }
	if(any(is.na(dcW.albc$countCells))){dcW.albc[which(is.na(dcW.albc$countCells)),]$countCells = 0 }
	dcW.albc=dcW.albc[-10,]
	dcW.albc$countCells=factor(dcW.albc$countCells, levels=c("< 10","10 - 20", "20 - 40", "40 - 60", ">60" ))
	dcW.albc$col=dcW.albc$countCells
	levels(dcW.albc$col) <- c("grey0", "grey14", "grey34", "grey64", "grey89")
	return(dcW.albc)
}

## barplot of the distribution of genetic diversity mean by cell with IC
plot_latbandGD <-function(dcW.lb, firstRow, secondRow,xlabFig,titleFig,labels.lb){	
	barplot.W <- barplot(dcW.lb[firstRow:secondRow,2],horiz=T,col=as.character(dcW.lb[firstRow:secondRow,5]),xlim=c(0,0.05),yaxs="i",cex.axis=1.2,tcl=-0.25,mgp=c(3,0.4,0),plot=TRUE,offset=0.000)
  	axis(side=2,at=barplot.W[firstRow:secondRow,1],labels=labels.lb,las=2,cex.axis=1.2,tcl=-0.25)
    arrows(dcW.lb[firstRow:secondRow,2]-dcW.lb[firstRow:secondRow,3],
    	barplot.W,
    	dcW.lb[firstRow:secondRow,2]+dcW.lb[firstRow:secondRow,3],
    	barplot.W,
    	angle=90, code=3,length=0.05)
	mtext(titleFig,side=2,line=0.4,at=21.5,las=2,cex=1.2)
	mtext(xlabFig,side=1,line=1.6,at=0.02,cex=0.85)
	box()
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

## pictures for decorum
Mraster_imag <- readPNG("01-infos/images/reef_fish.png")
Fraster_imag <- readPNG("01-infos/images/freshwater_fish.png")


##########################################################################
## add genetic diversity to the grid cells
gridM=grid_GD(grid_Wgs84,dcM.1,dcM)
gridF=grid_GD(grid_Wgs84,dcF.1,dcF)

## dataframe of genetic diversity by latband
dcM.lb=latbandGD(dcM.1,dcM,grid)
dcF.lb=latbandGD(dcF.1,dcF,grid)


##########################################################################
## write pdf files
labels.lb=c("-80","","-60","","-40","","-20","","0","","20","","40","","60","","80")

tiff(filename="10-figures/figure1.tiff",width=20,height=16,units="cm",res=640,compression="lzw")
#cairo_ps(filename="10-figures/figure1.eps",width=7.87,height=6.3,fallback_resolution=600,onefile=FALSE)

layout(mat=rbind(c(1,1,2),c(1,1,2),c(1,1,2),c(3,3,4),c(3,3,4),c(3,3,4)))
par(mar=c(3,3,3,1))
### marin
plot_GD_world(gridM,worldcoast,biglakes,riverlakes,"(a)")
plot_latbandGD(dcM.lb, 1,17,labels=labels.lb,"", "(b)")
rasterImage(Mraster_imag,xleft=0.038,ybottom=17.5,xright=0.049,ytop=20)

### freshwater
plot_GD_world(gridF,worldcoast,biglakes,riverlakes,"(c)")
legend("bottomleft",legend=levels(gridF$GD),col=levels(gridF$col),
       pch=15,cex=0.9,title="Genetic diversity",pt.cex=1)
plot_latbandGD(dcF.lb, 1,17,labels=labels.lb,"Genetic diversity", "(d)")
legend("bottomright",legend=levels(dcF.lb$countCells),col=levels(dcF.lb$col),
       pch=15,cex=0.9,title="# cells",pt.cex=1)
rasterImage(Fraster_imag,xleft=0.031,ybottom=17.5,xright=0.049,ytop=20)

dev.off()





