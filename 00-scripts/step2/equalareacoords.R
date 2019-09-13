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
############################################################################
## Get grid equal area ID in which a point (x,y) is located
##
## attributes at each sequence an ID of cell of the shapefile 
## of worldmap equal area projection from its coordinates.
##
############################################################################
library(rgeos)
library(rgdal)
## list of the cluster of species ".coords" files on which we want to get "equal area ID" of samples geo-position
IN_FILES = c("./06-species_alnt_cluster/total","./06-species_alnt_cluster/freshwater","./06-species_alnt_cluster/marine")
species_to_remove=c()
## convert ".coords files" of a species into "grid-equalarea ID" files
getEqualAreaSpecies <- function(nom_f){
	test = read.table(nom_f,header=F)	
	colnames(test)=c("lat","lon")
	test2 = cbind(test$lon,test$lat)
	testS = SpatialPoints(test2)
	proj4string(testS) <- CRS('+init=epsg:4326')
	SS =  spTransform(testS, proj4string(grid))
	coorS = SS@coords
	is_into=apply(coorS,1,function(x) grid@data$IDcell[which(grid@data$right > x[1] & grid@data$left < x[1] & grid@data$bottom < x[2] & grid@data$top > x[2])])
        if(length(is_into)!=0) {
	    #nouvelle stat coords
	    final_coords= na.omit(cbind(test$lon,test$lat,as.integer(is_into)))
	    new_nom_f=gsub("coords","equalareacoords",nom_f)
	    write.table(final_coords,new_nom_f,col.names=F,row.names=F,sep="\t")
        }
        else{
            print(nom_f)
            file.remove(gsub("coords","fasta",nom_f))
            file.remove(nom_f) 
        }
}

## apply getEqualAreaSpecies to each species in a cluster of species
equalAreaCluster <- function(liste_in_files){
    coords_files = Sys.glob(paste(liste_in_files,"/*.coords",sep=""))
    lapply(coords_files, function(x) getEqualAreaSpecies(x))
}

### load shapefile grid
grid <- readOGR("./01-infos/grid_equalarea200km","gridFish.b_260418")
grid$Id <- 1:nrow(grid@data)
### apply getEqualAreaSpecies to each cluster of species
lapply(IN_FILES, function (x) equalAreaCluster(x))
