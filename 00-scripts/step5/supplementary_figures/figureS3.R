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
## Supplementary Figure 3.
## Latitudinal distribution of species diversity
## for (a) boxplot of species diversity
## and (b) marine and (c) freshwater fish species
##
##########################################################################
## Libraries
lib_vect <-c("raster","plotrix","rgeos","rgdal","sp","maptools","shape","parallel","png","plyr")
sapply(lib_vect,library,character.only=TRUE)


##########################################################################
## functions
source("00-scripts/step5/figures/functions.R")

## convert coordinates into epsg:4326 and aggregate by latband(10degree)
## and gives the mean of the species diversity per species per cell per latband
latbandrichness <- function(dcW,gridW) {
    gridW.1=gridW[which(gridW$IDcell %in% dcW$cell),]	
	centroW <- gCentroid(gridW.1,byid=TRUE)	
	spW=SpatialPoints(data.frame(as.numeric(centroW@coords[,1]),as.numeric(centroW@coords[,2])))
	proj4string(spW)=proj4string(gridW)
	spW.g=spTransform(spW, CRS("+init=epsg:4326"))
	lat_W=matrix(spW.g@coords[,2])
	latband_W=apply(lat_W, 1, FUN=function(x) coords_to_latband(x))
	dcW.latband=data.frame(latband=latband_W,richness=dcW$richness)
	dcW.lb.mean=aggregate(richness ~ latband, FUN = mean, dcW.latband)
    names(dcW.lb.mean)=c("latband","richness_mean")
	latband=data.frame(latband=seq(from=-90,to=90,by=10)[-10])
    dcW.lb.sd=aggregate(richness ~ latband, FUN = function(x) sd(x)/sqrt(length(x)),dcW.latband)
    names(dcW.lb.sd)=c("latband","richness_sd")
	dcW.mnsd=merge(dcW.lb.mean,dcW.lb.sd,by="latband",all=T)
    dcW.alb=merge(latband, dcW.mnsd , by="latband",all=T)
	dcW.alb[is.na(dcW.alb$richness_mean),]$richness_mean = 0
	dcW.alb[is.na(dcW.alb$richness_sd),]$richness_sd = 0 
	#if(any(is.na(dcW.lb.sd$richness_sd))){dcW.lb.sd[which(is.na(dcW.lb.sd$richness_sd)),]$richness_sd = 0 }
	return(dcW.alb)
}

plot_latband_richness <-function(dcW.lb, firstRow, secondRow,xlabFig,titleFig){
	barplot.W <- barplot(dcW.lb[firstRow:secondRow,2],horiz=T,col="grey",xlim=c(0,1600),cex.axis=1.2,tcl=-0.25,mgp=c(3,0.4,0),plot=TRUE,offset=0.000,xaxt="n")
    axis(side=1,line=0,cex.axis=1.6,lwd=0.35,tcl=-0.25,bg="white",labels=as.character(seq(0,1600,by=400)), at=seq(0,1600,by=400),mgp=c(3,0.5,0))
	axis(side=2,at=barplot.W[firstRow:secondRow,1],labels=seq(-90,90,10)[-10],las=2,cex.axis=1.2,tcl=-0.25)
	arrows(dcW.lb$richness_mean-dcW.lb$richness_sd,
    	barplot.W,
    	dcW.lb$richness_mean+dcW.lb$richness_sd,
    	barplot.W,
    	angle=90, code=3,length=0.05)
    mtext(titleFig,side=2,line=1,at=23.5,las=2,cex=1.4)
	mtext(xlabFig,side=1,line=2,at=800,cex=1)
	box()
}

boxplot_richness_env <- function(gdW,titleFig,ylabFig) {
	boxplot(richness~environment,data=gdW, main="", xlab="", ylab="",col="grey",ylim=c(0,max(gdW$richness)),cex.axis=1.6,cex.lab=1,yaxt = "n")
	axis(side=2,line=0,cex.axis=1.6,lwd=0.35,tcl=-0.25,bg="white",labels=as.character(seq(0,3600,by=400)), at=seq(0,3600,by=400),mgp=c(3,0.5,0))
	mtext(titleFig,side=2,line=1,at=3850,las=2,cex=1.4)
	mtext(ylabFig,side=2,line=3,at=1600,cex=1.2)

}


##########################################################################
## load data
dc = load_data("09-all_descripteurs/total_data_genetic_diversity_with_all_descripteurs.tsv")
dcT=dc[[1]];dcM=dc[[2]];dcF=dc[[3]];dcM.n1=dc[[4]];dcF.n1=dc[[5]];dcM.1=dc[[6]];dcF.1=dc[[7]]
richness=read.csv("01-infos/equalarea_id_coordsCA_FWRS_MR_RS.csv")
## freshwater
Fr=merge(dcF.1,richness,by ="cell", all.x=F)
Fr=Fr[which(Fr$RS_FW_FISH !=0),]
## remove 0 value GD
Fr= Fr[which(Fr$GD_mean != 0),]
## filter wrong freshwater assignation
wrongFreshwaterCells=c(3304,4741,8453,9146,11400,11474,11548,12434,12730,11860)
Fr=Fr[which(!Fr$cell %in% wrongFreshwaterCells),]
Fr=Fr[,c(1,28)]
Fr=na.omit(Fr)
names(Fr)=c("cell","richness")
## marine
Mr=merge(dcM.1,richness,by ="cell", all.x=F)
Mr=Mr[which(Mr$RS_MR_FISH !=0),]
## remove 0 value GD
Mr= Mr[which(Mr$GD_mean != 0),]
Mr=Mr[,c(1,29)]
Mr=na.omit(Mr)
names(Mr)=c("cell","richness")
## total
Tr=rbind(Mr,Fr)
environment=c(rep("marine",dim(Mr)[1]),rep("freshwater",dim(Fr)[1]))
Tr=cbind(Tr,environment)
Tr$environment <- factor(Tr$environment,levels = c('marine','freshwater'),ordered = TRUE)


##########################################################################
## load worldcoast shapefile
grid <- readOGR("01-infos/grid_equalarea200km","gridFish.b_260418")


##########################################################################
## dataframe of genetic diversity by latband
dcM.lb=latbandrichness(Mr,grid)
dcF.lb=latbandrichness(Fr,grid)


##########################################################################
## write pdf files
pdf("10-figures/figureS3.pdf",width=17,height=8,paper='special')

layout(matrix(c(1,2,3),nrow=1))
par(mar=c(4,5,4,4))


boxplot_richness_env(Tr,"(a)","Species diversity")
plot_latband_richness(dcM.lb, 1,18,"Species diversity", "(b)")
plot_latband_richness(dcF.lb, 1,18,"Species diversity", "(c)")

dev.off()
