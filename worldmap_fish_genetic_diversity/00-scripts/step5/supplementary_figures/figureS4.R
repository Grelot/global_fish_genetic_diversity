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
## Figure S4 : Regional effect on the global genetic diversity pattern
## (a) for marine and  (b) freshwater fish species.
##
##
##########################################################################
## Libraries
lib_vect <-c("MASS","hier.part","countrycode","lme4","sjPlot","raster","plotrix","rgeos","rgdal","sp","maptools","shape","parallel","png","plyr")
sapply(lib_vect,library,character.only=TRUE)


##########################################################################
## functions
source("00-scripts/step5/figures/functions.R")

boxplot_gd_region <- function(gdW,titleFig,ylabFig) {
	boxplot(gd~regions,data=gdW, main="", xlab="", ylab="",col="grey",ylim=c(0,0.1),cex.axis=1.4,cex.lab=1.4)
	mtext(titleFig,side=2,line=0.4,at=4.6,las=2,cex=1.6)
	mtext(ylabFig,side=2,line=2.5,at=0.05)
}


##########################################################################
## load data
dc = load_data("09-all_descripteurs/total_data_genetic_diversity_with_all_descripteurs.tsv")
dcT=dc[[1]];dcM=dc[[2]];dcF=dc[[3]];dcM.n1=dc[[4]];dcF.n1=dc[[5]];dcM.1=dc[[6]];dcF.1=dc[[7]]


##########################################################################
## MARINE
### add region indopacific/atlantic information for MARINE 
lon = as.double(dcM.n1$x)
lat = as.double(dcM.n1$y)
lonlat=cbind(lon,lat)
lonlat_s=SpatialPoints(lonlat)
atlantic_y=c(34.5,51.5,80,85,85,85,36.5,27,16,14.25,11.68,8.8,9.17,6.91,-4.4,-15.5,-17.7,-26.1,-51,-55,-64.1,-82,-82,20,30)
atlantic_x=c(39.5,49.3,115,21,77,-119,-108,-105,-95,-87.5,-84.7,-82.4,-79,-75.7,-74.6,-70.4,-68.8,-68.2,-72.2,-66,-60,-52.5,21,21,33.5)
atlantic_xym=cbind(atlantic_x,atlantic_y)
atlantic_p = Polygon(atlantic_xym)
atlantic_ps = Polygons(list(atlantic_p),1)
atlantic_sps = SpatialPolygons(list(atlantic_ps))
grid <- readOGR("./01-infos/grid_equalarea200km","gridFish.b_260418")
proj4string(atlantic_sps) <- CRS("+proj=longlat")
proj4string(lonlat_s) <- CRS(proj4string(grid))
atlantic_sps = spTransform(atlantic_sps, proj4string(grid))
atlantic_id=which(over(lonlat_s,atlantic_sps)==1)
indopacific_id=which(is.na(over(lonlat_s,atlantic_sps)))
dcM.atlantic=dcM.n1[atlantic_id,]
dcM.indopac=dcM.n1[indopacific_id,]
regions=c(rep("Atlantic",length(dcM.atlantic$y)),rep("Indopacific",length(dcM.indopac$y)))
dcM.r=rbind(dcM.atlantic,dcM.indopac)
dcM.r=cbind(dcM.r,regions)

### format data
dcM.r1=format_dc(dcM.r)
names(dcM.r1)

dcM.r2=data.frame(
	cell=dcM.r1$cell,
	rID=c(dcM.r1$regions),
	regions=dcM.r1$regions,
	GD_mean=dcM.n1[which(dcM.n1$cell %in% dcM.r1$cell),]$GD_mean)

dcM.r3<-na.omit(dcM.r2)

dM_GD.r=data.frame(gd=dcM.r3$GD_mean,regions=dcM.r3$regions)


##########################################################################
## FRESHWATER
### add region continent information for FRESHWATER
north_america=c("USA","MEX","CAN")
regions=countrycode(dcF.1$ISO3,"iso3c","continent")
regions[which(dcF.1$ISO3 %in% north_america)]="North america"
regions[which(regions=="Americas")]="South america"
regions[which(is.na(regions))]="Antarctica"
dcF.r1=cbind(dcF.1,regions)

### format data
dcF.r2=data.frame(
	cell=dcF.r1$cell,
	rID=c(dcF.r1$regions),
	regions=dcF.r1$regions,
	GD_mean=dcF.n1[which(dcF.n1$cell %in% dcF.r1$cell),]$GD_mean)
dcF.r3=na.omit(dcF.r2)
#Additional filters : Remove oceania and Antartica and bad anotated cells and run model again : 346 cells
dcF.r3=dcF.r3[dcF.r3$rID != 6,]
dcF.r3=dcF.r3[dcF.r3$rID != 2,]
dcF.r3=dcF.r3[dcF.r3$cell!=3304,]
dcF.r3=dcF.r3[dcF.r3$cell!=9146,]

dF_GD.r=data.frame(gd=dcF.r3$GD_mean,regions=droplevels(dcF.r3$regions))
droplevels(dF_GD.r$regions)

##########################################################################
## write pdf files
pdf("10-figures/figureS4.pdf",width=14,height=8,paper='special')

layout(matrix(c(1,1,2,2,2),nrow=1))
par(mar=c(4,4,4,2))
boxplot_gd_region(dM_GD.r,"(a)","Genetic diversity")
boxplot_gd_region(dF_GD.r,"(b)","")

dev.off()
