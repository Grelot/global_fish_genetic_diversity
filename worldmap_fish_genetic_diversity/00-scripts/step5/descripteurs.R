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
## Generates a table with row as cell ID into a shapefile grid nested equal area
## Each column is a descriptor :
## - Marine oxygen concentration
## - Sea Surface Temperature (marine) or Global Mean Temperature (freshwater)
## - Velocity
## - Drainage basin surface area
##########################################################################

## Libraries
library(raster)
library(sp)
library(rgdal)
library(rgeos)

## Functions
# coordonnate x,y of each cell
center_of_square <- function(tempo) {
    y=(tempo$bottom+tempo$top)/2
    x=(tempo$right+tempo$left)/2
    tempo = tempo[,-which(names(tempo) %in% c("left","right","top","bottom"))]
    tempo= cbind(x,y,tempo)
    colnames(tempo)[1] = "x"
    colnames(tempo)[2]="y"
    return(tempo)
}

# from O2, SST, VEL,basin, RASTERS and metrics by area FILE, generate a dataframe that can be used for lm
dataframe_metrics <- function(grid,basin,O2,SST,VEL,metrics_by_area, buffer_distance) {
    cells = center_of_square(metrics_by_area)
    xy = data.frame(cbind(cells$x, cells$y))
    colnames(xy) = c("x","y")
    coordinates(xy) = c("x", "y")
    proj4string(xy) = CRS(proj4string(grid))
    xy = spTransform(xy,CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
    basin.proj=spTransform(basin,proj4string(xy))
    basin.extr=extract(basin.proj,xy,buffer=buffer_distance,fun=min)
    Env = data.frame(
          #O2=extract(O2, xy,buffer=buffer_distance,fun=mean),
          #SST=extract(SST, xy,buffer=buffer_distance,fun=mean),
          #VEL=extract(VEL, xy,buffer=buffer_distance,fun=mean),
          basin_area=basin.extr$Surf_area)
    Final = data.frame(cbind(cells, Env))
    rownames(Final) <- NULL
    return(Final)
}

## Spatial shapefile
grid <- readOGR("./01-infos/grid_equalarea200km","gridFish.b_260418")
basin <- readOGR("./01-infos/datatoFigshare/","Basin042017_3119")


## Spatial layers
# MARINE Oxygen concentration [?mol/l]
marine_O2 = raster("01-infos/spatial_layers/marine_bo_o2dis.asc",crs="+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
# MARINE Sea surface temperature [C?]
marine_SST = raster("01-infos/spatial_layers/marine_bo_sst_mean.asc",crs="+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
# MARINE Velocity
marine_VEL = raster("01-infos/spatial_layers/marine_velocity_mean.asc",crs="+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
# FRESHWATER Global Mean Temperature
freshwater_GMT=raster("01-infos/spatial_layers/freshwater_wc2.0_bio_10m_01.tif")
# FRESHWATER Velocity (behrman projection)
freshwater_VEL=raster("01-infos/spatial_layers/freshwater_velocity_mean.tif")

## GloRiC layers
#gloric=readOGR("./01-infos/GloRiC_v10_shapefile/","GloRiC_v10")
#saveRDS(gloric, "09-all_descripteurs/gloric.rds")
#gloric = readRDS("09-all_descripteurs/gloric.rds")

# buffer distance
buffer_distance=504000

##charger les donnees
# descripteurs environmentaux
dataCell=read.table("01-infos/datacell_grid_descriteurs.csv",header=T,sep=",")
# genetic diversity (marine)
metricsMarine=read.table("08-genetic_diversity/metrics_by_area_marine.csv",header=T,sep=",")
dcMm=merge(metricsMarine, dataCell, by="cell")
# descripteurs climatiques (marine)
dcM=dataframe_metrics(grid,basin,marine_O2,marine_SST,marine_VEL,dcMm ,buffer_distance)
# genetic diversity (freshwater)
metricsFreshwater=read.table("08-genetic_diversity/metrics_by_area_freshwater.csv",header=T,sep=",")
dcFm=merge(metricsFreshwater, dataCell, by="cell")
# descripteurs climatiques (freshwater)
dcF=dataframe_metrics(grid,basin,marine_O2,freshwater_GMT,freshwater_VEL,dcFm,buffer_distance)
# clustering marine/freshwater
is_marine=c(rep(0,dim(dcF)[1]),rep(1,dim(dcM)[1]))
dcT=rbind(dcF,dcM)
dcTc=cbind(dcT,is_marine)

##write file with all data cell genetic diversity with all descripteurs
write.table(dcTc,"09-all_descripteurs/total_data_genetic_diversity_with_all_descripteurs.tsv",col.names=T,row.names=F,sep="\t")


