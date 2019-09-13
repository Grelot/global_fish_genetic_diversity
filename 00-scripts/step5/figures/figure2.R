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
##
## Figure 2 : Congruence between fish genetic and species diversity.
## Outputs of the linear model (lm) testing the correlation between genetic diversity and marine
## (a) and freshwater (c) species diversity.
## Values of genetic diversity were reported on the global map using a color
## gradient depending on the values of the correlation 
## for marine species (b) and freshwater species (d) respectively
##
##
##
##########################################################################
## Libraries
lib_vect <-c("plyr","raster","sp","rgeos","rgdal","sp","ggplot2","MASS","SpatialPack","maptools")
sapply(lib_vect,library,character.only=TRUE)


##########################################################################
## functions
source("00-scripts/step5/figures/functions.R")

## from (x,y) values, it generates a spatial color
fun_xy <- function(x, y) {
  R <- (x+1)/2
  G <- (1-x)/2
  B <- (y+1)/2
  A <- 1- 0.5*exp(-(x^2+y^2)/0.2)
  rgb(R, G, B, A)
}

## rotate a matrix "x" in clockwise order
rotate_clockwise <- function(x) { t(apply(x, 2, rev))}


## COULEUR POUR LES DALTONIENS
## col <- rev(c("#d73027","#f46d43","#fdae61","#fee090","#ffffbf","#e0f3f8","#abd9e9","#74add1","#4575b4"))

## return an integer (1,2,3,4) corresponding to the quantile of the val
quantile_value <- function(val) {
	valQuantile=quantile(val)
	valIntervals=as.numeric(cut(val,valQuantile))
	valIntervals[is.na(valIntervals)] <- 1
	return(valIntervals)
}

## from GD and species diversity value, attribute a color in a table 4*4
## color is assigned according to the quantile of the value [GD,richness]
color_value <- function(GD, richness) {
	#matCol = outer(seq(-1,1,length=4), seq(-1,1,length=4), FUN = fun_xy)
    matCol=matrix( 
    c("#fee090","#fdae61","#f46d43","#d73027",
        "#F4E6B3","#E2BC8E","#C98272","#A64756",
        "#EAEDD5","#C6CBBC","#9F98A2","#765E85",
        "#e0f3f8","#abd9e9","#74add1","#4575b4"),
        nrow=4,              # number of rows 
        ncol=4,              # number of columns 
        byrow = FALSE)
	matCol= rotate_clockwise(t(matCol))
	GDInterval=quantile_value(GD)
	richnessInterval=quantile_value(richness)
	colVal=rep("NA",length(GDInterval))
	for(i in seq(1,length(colVal))) {
		colVal[i]=matCol[GDInterval[i],richnessInterval[i]]
	}
	return(colVal)
}

## build figure from a linear model of lines and points with ggplot
## it requires title and x lab and y lab arguments
## it requires linear model and x and y data
## assign a 2D-colour to each point according to their value (x,y)
plot_color_lm<- function(mod,df_x,df_y,x_lab,y_lab,titleFig) {
    #R square
    rqs = summary(mod)$adj.r.squared
    # correlation
    corrValue=round(sqrt(rqs),4)
    #coefficients
    b = summary(mod)$coefficients[1]
    a1 = summary(mod)$coefficients[2]    
    df=data.frame(df_x,df_y)    
    dr_x=df_x
    dr_y=b+a1*dr_x
    dregr=data.frame(dr_x,dr_y)
    plot(df_x,df_y,col=color_value(df_x,df_y),pch=19,cex=0.65,lwd=2,axes=FALSE,xlab="",ylab="",xlim=c(0,0.065),ylim=c(0,3600))
    lines(dr_x,dr_y,lwd=2)
    axis(side=1,line=0,cex.axis=0.85,lwd=0.35,tcl=-0.25,bg="white",labels=as.character(seq(0,0.6,by=0.02)), at=seq(0,0.6,by=0.02),mgp=c(3,0.5,0))
    axis(side=2,line=0,cex.axis=0.85,lwd=0.35,tcl=-0.25,bg="white",labels=as.character(seq(0,3600,by=400)), at=seq(0,3600,by=400),mgp=c(3,0.5,0))
    box()
    mtext(x_lab,side=1,line=1.5,cex=0.7)
    mtext(y_lab,side=2,line=1.5,cex=0.7)
    mtext(paste("r =",corrValue),side=3,line = -1.5,at=0.05,cex=0.7)
    mtext(titleFig,line=0.4,at=0)
}

## add GD and species diversity information and color to the grid
grid_GD_richness <- function(gridW,dcW) {	
	gridW.1=gridW[which(gridW$IDcell %in% dcW$cell),]	
	centroW <- gCentroid(gridW.1,byid=TRUE)
	dfW=data.frame(lon=as.numeric(centroW@coords[,1]), lat=as.numeric(centroW@coords[,2]),col=color_value(dcW$GD,dcW$richness))
	return(dfW)
}

## world map with colored squares according to their (richness,genetic div) values
plot_GD_world <- function(gridW,worldcoast,biglakes,riverlakes,titleFig) {
    plot_world(worldcoast,biglakes,riverlakes)
    points(gridW$lon,gridW$lat,bg=as.character(gridW$col),pch=22,cex=0.65,col="black",lwd=0.2)
    box()
    mtext(titleFig,side=2,line=0.4,at=90,las=2)
}

## plot a legend color 2D
color_map_legend <- function (colours, labels = FALSE, borders = NULL, cex_label = 1,xlab= NULL,ylab=NULL) {
    n <- length(colours)
    ncol <- ceiling(sqrt(n))
    nrow <- ceiling(n/ncol)
    colours <- c(colours, rep(NA, nrow * ncol - length(colours)))
    colours <- matrix(colours, ncol = ncol, byrow = TRUE)
    old <- par(pty = "s", mar = c(0, 0, 0, 0))
    on.exit(par(old))
    size <- max(dim(colours))
    plot(c(0, size), c(0, -size), type = "n",axes=F)
    rect(col(colours) - 1, -row(colours) + 1, col(colours), -row(colours), 
        col = colours, border = borders)
    if (labels) {
        text(col(colours) - 0.5, -row(colours) + 0.5, colours, 
            cex = cex_label)
    }
    mtext(xlab,side=1,line=0.2,cex=0.4)
    mtext(ylab,side=2,line=0.2,cex=0.4)
}


##########################################################################
## load data
dc = load_data("09-all_descripteurs/total_data_genetic_diversity_with_all_descripteurs.tsv")
dcT=dc[[1]];dcM=dc[[2]];dcF=dc[[3]];dcM.n1=dc[[4]];dcF.n1=dc[[5]];dcM.1=dc[[6]];dcF.1=dc[[7]]
richness=read.csv("01-infos/equalarea_id_coordsCA_FWRS_MR_RS.csv")
## freshwater
Fr=merge(dcF.n1,richness,by ="cell", all.x=F)
Fr=Fr[which(Fr$RS_FW_FISH !=0),]
## remove 0 value GD
Fr= Fr[which(Fr$GD_mean != 0),]
## filter wrong freshwater assignation
wrongFreshwaterCells=c(3304,4741,8453,9146,11400,11474,11548,12434,12730,11860)
Fr=Fr[which(!Fr$cell %in% wrongFreshwaterCells),]
Fr=Fr[,c(1,2,3,4,6,7,28)]
Fr=na.omit(Fr)
names(Fr)=c("cell","x","y","GD_mean","nb_species","nb_indv_mean","richness")
## marine
Mr=merge(dcM.n1,richness,by ="cell", all.x=F)
Mr=Mr[which(Mr$RS_MR_FISH !=0),]
Mr=Mr[,c(1,2,3,4,6,7,29)]
Mr=na.omit(Mr)
Mr= Mr[which(Mr$GD_mean != 0),]
names(Mr)=c("cell","x","y","GD_mean","nb_species","nb_indv_mean","richness")


##########################################################################
## load worldcoast shapefile
grid <- readOGR("01-infos/grid_equalarea200km","gridFish.b_260418")
worldcoast <- readOGR("01-infos/ne_50m_land",layer="ne_50m_land")
## load rivers and lake centerlines shapefile
riverlakes  <- readOGR("01-infos/ne_50m_rivers_lake_centerlines_scale_rank",layer="ne_50m_rivers_lake_centerlines_scale_rank")
grid_Wgs84 <- spTransform(grid, proj4string(worldcoast))
biglakes <- readShapePoly("01-infos/big_lakes/GSHHS_h_L2.shp",proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"))


##########################################################################
## models & statistical tests
## linear regression ~ correlation
### step AIC
modFr=lm(GD_mean~richness+nb_species+nb_indv_mean,data=Fr)
stepAIC(modFr)
### final model
modFr.final=lm(richness~GD_mean,data=Fr)
#modFr.t=lm(scale(log(GD_mean))~scale(log(richness)),data=Fr)
modMr.final=lm(richness~GD_mean,data=Mr)
##### Modified t-test of spatial association test
###### marine
coordsM = cbind(Mr$x,Mr$y)
zM <- cor.spatial(Mr$GD_mean,Mr$richness,coordsM)
modified.ttest(Mr$GD_mean,Mr$richness,coordsM, nclass = 13) 

###### freshwater
coordsF = cbind(Fr$x,Fr$y)
modified.ttest(Fr$GD_mean,Fr$richness,coordsF, nclass = 13) 


##########################################################################
### grid with species diversity and genetic diversity values
grid.M=grid_GD_richness(grid_Wgs84,Mr)
grid.F=grid_GD_richness(grid_Wgs84,Fr)

## legend matrix of colors
#legMatCol = outer(seq(-1,1,length=4), seq(-1,1,length=4), FUN = fun_xy)
legMatCol=matrix( 
 c("#fee090","#fdae61","#f46d43","#d73027",
    "#F4E6B3","#E2BC8E","#C98272","#A64756",
    "#EAEDD5","#C6CBBC","#9F98A2","#765E85",
    "#e0f3f8","#abd9e9","#74add1","#4575b4"),
    nrow=4,              # number of rows 
    ncol=4,              # number of columns 
    byrow = FALSE)


##########################################################################
## write pdf files
#tiff(filename="10-figures/figure2_colorblind.tiff",width=20,height=16,units="cm",res=640,compression="lzw")
cairo_ps(filename="10-figures/figure2_colorblind.eps",width=7.87,height=6.3,fallback_resolution=600,onefile=FALSE)

layout(mat=rbind(c(1,2,2),c(1,2,2),c(1,2,2),c(3,4,4),c(3,4,4),c(3,4,4)))
### marin
par(mar=c(2.5,3.5,2,2))
plot_color_lm(modMr.final,Mr$GD_mean,Mr$richness,"Genetic diversity","Species diversity","(a)")
par(mar=c(2.5,2,2,2))
plot_GD_world(grid.M,worldcoast,biglakes,riverlakes,"(b)")
### freshwater
par(mar=c(2.5,3.5,2,2))
plot_color_lm(modFr.final,Fr$GD_mean,Fr$richness,"Genetic diversity","Species diversity","(c)")
par(mar=c(2.5,2,2,2))
plot_GD_world(grid.F,worldcoast,biglakes,riverlakes,"(d)")
## legend
#par( fig=c(0.4,0.45,0.15,0.2),new = T)
par( fig=c(0.39,0.46,0.11,0.18),new = T)
color_map_legend(legMatCol,xlab="Genetic Diversity\nQuantiles",ylab="Species diversity\nQuantiles")

dev.off()
