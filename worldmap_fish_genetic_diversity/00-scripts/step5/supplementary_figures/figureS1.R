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
## Supplementary Figure S1: Spatial autocorrelogramme based on the I-Moran coefficient
## (R package pgirmess function correlog )of the genetic diversity for
## (a) freshwater species.
## (b) saltwater species
## Distances classes are in km.
## Red dots indicated statistically significant values (p<0.05).
##
##
##########################################################################
## Libraries
lib_vect <-c("pgirmess","ggplot2","raster","rgdal","rgeos","MASS","hier.part","countrycode","lme4","sjPlot","plyr","gridExtra","visreg")
sapply(lib_vect,library,character.only=TRUE)

##########################################################################
## functions
source("00-scripts/step5/figures/functions.R")

cell_coordinates_meters_kilometers <- function(dcW.n1,dcW.r3) {
## cell coordinates in meters and kilometers
	coordsWm=data.frame(lon=dcW.n1$x,lat=dcW.n1$y,cell=dcW.n1$cell)
    coordsWm.r3=coordsWm[which(coordsWm$cell %in% dcW.r3$cell),c(1,2)]
    coordsWkm.r3<-coordsWm.r3/1000
    return(coordsWkm.r3)
}

##########################################################################
## load data
dc = load_data("09-all_descripteurs/total_data_genetic_diversity_with_all_descripteurs.tsv")
dcT=dc[[1]];dcM=dc[[2]];dcF=dc[[3]];dcM.n1=dc[[4]];dcF.n1=dc[[5]];dcM.1=dc[[6]];dcF.1=dc[[7]]

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
regions=c(rep("atlantic",length(dcM.atlantic$y)),rep("indopacific",length(dcM.indopac$y)))
dcM.r=rbind(dcM.atlantic,dcM.indopac)
dcM.r=cbind(dcM.r,regions)

### format data
dcM.r1=format_dc(dcM.r)
dcM.r2<-dcM.r1[,-c(5,8,9,10,11,13, 14,16,17,18,21,22,23)]
dcM.r3<-na.omit(dcM.r2)

## cell coordinates in meters and kilometers
coordsMkm.r3=cell_coordinates_meters_kilometers(dcM.n1,dcM.r3)

## Spatial autocorrelogramme of the variable GD and for model residuals
#plot(coordsMkm.r3$lon,coordsMkm.r3$lat)
### Genetic diversity Spatial Autocorrelation Marin
#### coor center to 0
pgi.corM <- correlog(coords=coordsMkm.r3, z=dcM.r3$GD_mean, method="Moran", nbclass=100)
#plot(pgi.corM,xlim=c(0,10000))


##########################################################################
## FRESHWATER
### add region continent information for FRESHWATER
north_america=c("USA","MEX","CAN")
regions=countrycode(dcF.1$ISO3,"iso3c","continent")
#regions=countrycode(dcF$ISO3,"iso3c","continent")
regions[which(dcF.1$ISO3 %in% north_america)]="North_america"
#regions[which(dcF$ISO3 %in% north_america)]="North_america"

regions[which(regions=="Americas")]="South_america"
regions[which(is.na(regions))]="Antarctica"
dcF.r1=cbind(dcF.1,regions)

## format data
dcF.r2=data.frame(x=dcF.r1$x,
	y=dcF.r1$y,
	cell=dcF.r1$cell,
	regions=c(dcF.r1$regions),
	GD_mean=dcF.r1$GD_mean)
dcF.r3=na.omit(dcF.r2)
#Additional filters : Remove oceania and Antartica and bad anotated cells and run model again : 346 cells
dcF.r3=dcF.r3[dcF.r3$regions != 6,]
dcF.r3=dcF.r3[dcF.r3$regions != 2,]
dcF.r3=dcF.r3[dcF.r3$cell!=3304,]
dcF.r3=dcF.r3[dcF.r3$cell!=9146,]


## cell coordinates in meters and kilometers
coordsFkm.r3=cell_coordinates_meters_kilometers(dcF.n1,dcF.r3)

## Residual models spatial autocorrelation
pgi.corF <- correlog(coords=coordsFkm.r3, z=dcF.r3$GD_mean, method="Moran", nbclass=100)
#plot(pgi.corF,xlim=c(0,10000))


##########################################################################
## write pdf files
pdf("10-figures/figureS1.pdf",width=14,height=10,paper='special')

layout(matrix(c(1,1,2,2),nrow=1))
par(mar=c(4,4,4,2)) # c(bottom, left, top, right)
plot(pgi.corF,xlim=c(200,10000),ylim=c(-0.3,0.3),main="", xlab="", ylab="",cex.axis=1.6,xaxt="n")
axis(side=1, at=c(200,2000,4000,6000,8000,10000), labels=c(200,2000,4000,6000,8000,10000),mgp=c(3,1,0),cex.axis=1.6)

mtext("(a)",side=3,line=1,at=0,cex=1.6)
mtext("Moran I",side=2,line=2,at=0.05,cex=1.4)
mtext("Distance classes (km)",side=1,line=3,cex=1.4)

plot(pgi.corM,xlim=c(0,10000),ylim=c(-0.3,0.3),main="", xlab="", ylab="",cex.axis=1.6,xaxt="n")
axis(side=1, at=c(200,2000,4000,6000,8000,10000), labels=c(200,2000,4000,6000,8000,10000),mgp=c(3,1,0),cex.axis=1.6)
mtext("(b)",side=3,line=1,at=0,cex=1.6)
mtext("Distance classes (km)",side=1,line=3,cex=1.4)

dev.off()