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
## Supplementary Figure S8: 
##
##
## We calculated intraspecific genetic diversity mean in each 10° latitudinal bands.
## To estimate confidence interval of this mean intraspecific genetic diversity of each latitudinal band, we did 1000 boostrap replicates. 
## We resampled the intraspecific genetic diversity from each species and we calculated the standard deviation.
##
##
## (a) saltwater species
## (b) freshwater species.
##
## Genetic diversity by latband
##
##
##########################################################################
## Libraries

library("png")
##########################################################################
## load data


## from latbands_numbers.csv file return a dataframe
## latband GD_band_mean   GD_band_sd
latband_GD <- function(ltbdFile) {
    ltbdNumbers=read.table(ltbdFile,header=T,sep=",")
    names(ltbdNumbers)=c("latband","GD_band_mean","GD_band_sd")
    latband=data.frame(latband=seq(from=-90,to=90,by=10))
    dcW.alb=merge(latband, ltbdNumbers , by="latband",all=T)
    for(i in seq(1,length(latband$latband)-1)) {
        dcW.alb$latband[i]=paste(latband$latband[i],"° - ", latband$latband[i]+10,"°",sep="")
    }
    dcW.albf=data.frame(latband=dcW.alb$latband,
    GD_band_mean=dcW.alb$GD_band_mean,
    GD_band_sd=dcW.alb$GD_band_sd)
    dcW.albf=dcW.albf[-19,]
    return(dcW.albf)
}

## barplot of the distribution of genetic diversity mean by cell with IC
plot_latbandGD <-function(dcW.lb, firstRow, secondRow,xlabFig,titleFig){
	barplot.W <- barplot(dcW.lb[firstRow:secondRow,2],horiz=T,col="grey",xlim=c(0,0.05),cex.axis=0.95,tcl=-0.25,mgp=c(3,0.4,0),plot=TRUE,offset=0.000)
  	axis(side=2,at=barplot.W[1:18,1],labels=dcW.lb[1:18,1],las=2,cex.axis=0.95,tcl=-0.25)
    arrows(dcW.lb$GD_band_mean-dcW.lb$GD_band_sd,
    	barplot.W,
    	dcW.lb$GD_band_mean+dcW.lb$GD_band_sd,
    	barplot.W,
    	angle=90, code=3,length=0.05)
	mtext(titleFig,side=2,line=0.4,at=22.5,las=2)
	mtext(xlabFig,side=1,line=1.6,at=0.02,cex=0.75)
	box()
}

##########################################################################
## pictures for decorum
Mraster_imag <- readPNG("01-infos/images/reef_fish.png")
Fraster_imag <- readPNG("01-infos/images/freshwater_fish.png")

##########################################################################
## FRESHWATER
latbandsFile="08-genetic_diversity/freshwater_latbands_bootstraps.csv"
ltbdF=latband_GD(latbandsFile)
## filter below -40 latitude (we removed these species from the analysis)
ltbdF[c(2,3),c(2,3)]=NA


## MARINE
latbandsFile="08-genetic_diversity/marine_latbands_bootstraps.csv"
ltbdM=latband_GD(latbandsFile)

##########################################################################
## write pdf files
pdf("10-figures/figureS8.pdf",width=14,height=9,paper='special')

layout(matrix(c(1,1,2,2),nrow=1))
par(mar=c(4,6,4,2))
plot_latbandGD(ltbdM, 1,18, "Genetic Diversity", "(a)" )
rasterImage(Mraster_imag,xleft=0.032, ybottom=19.5,xright =0.042 ,ytop=22)
plot_latbandGD(ltbdF, 1,18, "Genetic Diversity", "(b)" )
rasterImage(Fraster_imag,xleft=0.023, ybottom=19.5,xright =0.042 ,ytop=22)


dev.off()