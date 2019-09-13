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
## Supplementary Figure 7. Taxonomic coverage estimated in each cell
## as the (number of species used to estimate genetic diversity)/species richness
## for (a) marine species (b) freshwater species.
##
##########################################################################
## Libraries
lib_vect <-c("raster","plotrix","rgeos","rgdal","sp","maptools","shape","parallel","png","plyr")
sapply(lib_vect,library,character.only=TRUE)


##########################################################################
## functions
source("00-scripts/step5/figures/functions.R")

## convert number of sequences value into color value
color_taxcov_spe <- function(gd){
	gdc=gd
	gdc[which(gd<0.01)] = "0 - 1"	
	gdc[which(gd>=0.01 & gd<0.02)] = "1 - 2"
	gdc[which(gd>=0.02 & gd<0.05)] ="2 - 5"	
	gdc[which(gd>=0.05)] =">5"
	return(gdc)
}

## from a grid and table of descriptors
## create a grid with number of sequences information for each cell of the grid
grid_numberofseq <- function(gridW,dcW.1,dcW) {
    dcW.0=dcW[which(dcW$cell %in% dcW.1$cell),]
	specW= dcW.0$nb_species
    taxonSpe=dcW.1[which(dcW.1$cell %in% dcW.0$cell),]$richness
    taxCoverageSpe=specW/taxonSpe
	gridW.1=gridW[which(gridW$IDcell %in% dcW.1$cell),]	
	centroW <- gCentroid(gridW.1,byid=TRUE)
	dfW=data.frame(lon=as.numeric(centroW@coords[,1]),
	 lat=as.numeric(centroW@coords[,2]),
     taxonCoverageSpecies=color_taxcov_spe(taxCoverageSpe),
     taxonCoverageSpecies_color=color_taxcov_spe(taxCoverageSpe),
     tc_spe=taxCoverageSpe)
    dfW.ordered=dfW[order(dfW$tc_spe),]
	dfW.ordered$taxonCoverageSpecies=factor(dfW.ordered$taxonCoverageSpecies,
	levels=c("0 - 1","1 - 2", "2 - 5", ">5" ))
	dfW.ordered$taxonCoverageSpecies_color=dfW.ordered$taxonCoverageSpecies
    dfW.ordered=dfW.ordered[order(dfW.ordered$tc_spe),]
	levels(dfW.ordered$taxonCoverageSpecies_color) <- c("#FFCCCC","#FF6666","#FF0000","#990000")
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
## total number of species by cell
richness=read.csv("01-infos/equalarea_id_coordsCA_FWRS_MR_RS.csv")
## freshwater
Fr=merge(dcF.1,richness,by ="cell", all.x=F)
Fr=Fr[which(Fr$RS_FW_FISH !=0),]
## remove 0 value GD
Fr= Fr[which(Fr$GD_mean != 0),]
Fr=Fr[,c(1,2,3,4,6,7,28)]
Fr=na.omit(Fr)
names(Fr)=c("cell","x","y","GD_mean","nb_species","nb_indv_mean","richness")
## marine
Mr=merge(dcM.1,richness,by ="cell", all.x=F)
Mr=Mr[which(Mr$RS_MR_FISH !=0),]
Mr=Mr[,c(1,2,3,4,6,7,29)]
Mr=na.omit(Mr)
Mr= Mr[which(Mr$GD_mean != 0),]
names(Mr)=c("cell","x","y","GD_mean","nb_species","nb_indv_mean","richness")


##########################################################################
## taxon coverage : number of cells above a threshold of taxon coverage (level: species)
numberCellSpeTaxonCovM=c(0,0,0,0)
numberCellSpeTaxonCovF=c(0,0,0,0)
i=1
for(t_tc in c(0,0.01,0.02,0.05)) {
    numberCellSpeTaxonCovM[i]=length(which(Mr$nb_species/Mr$richness >=t_tc))
    numberCellSpeTaxonCovF[i]=length(which(Fr$nb_species/Fr$richness >=t_tc))
    i=i+1
}
numberCellSpeTaxonCov=rbind(numberCellSpeTaxonCovM,numberCellSpeTaxonCovF)
colnames(numberCellSpeTaxonCov)= c("no filter", "> 0.01",">0.02",">0.05")
row.names(numberCellSpeTaxonCov)=c("marine","freshwater")
##########################################################################
## write CSV files
write.csv(numberCellSpeTaxonCov,file="12-taxonomic_coverage/species_taxon_coverage_number_cells.csv",sep=",")



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
gridM.1.df=grid_numberofseq(grid_Wgs84,Mr,dcM)
gridF.1.df=grid_numberofseq(grid_Wgs84,Fr,dcF)

##########################################################################
## write pdf files
tiff(filename="10-figures/figureS7.tiff",width=26,height=22,units="cm",res=640,compression="lzw")
layout(mat=rbind(1,2))
par(mar=c(1.5,2,2,1))

## taxonomic coverage species cell
plot_nseq_world(gridM.1.df,4,worldcoast,biglakes,riverlakes,"(a)")
plot_nseq_world(gridF.1.df,4,worldcoast,biglakes,riverlakes,"(b)")
legend("bottomleft",legend=levels(gridM.1.df$taxonCoverageSpecies),col=levels(gridM.1.df$taxonCoverageSpecies_color),
       pch=15,cex=0.9,title="% taxon cover",pt.cex=1)

dev.off()
