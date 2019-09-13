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
## 
## Sensitivity analysis. impact of the amount of genetic data on model outputs.
## keep only grid cells with top values.
## filtering cells by number of individuals and number of available species per cell
## 
##
##
##########################################################################
## Libraries
lib_vect <-c("SpatialPack","car","ggplot2","raster","rgdal","rgeos","MASS","hier.part","countrycode","lme4","sjPlot","plyr","gridExtra","visreg")
sapply(lib_vect,library,character.only=TRUE)

##########################################################################
## functions
source("00-scripts/step5/figures/functions.R")

## autocorrelative variable
autocorr_var <- function(Wr,dcW,gridW,neighbourhoodRadius){
    GD_mean=c(Wr$GD_mean)
    coordNS= data.frame(cbind(dcW$x,dcW$y,dcW$cell))
    coordNS=coordNS[which(coordNS$X3 %in% Wr$cell),]
    coordsW<-data.frame(cbind(coordNS$X1, coordNS$X2)) #Coordinates in meter
    #coordsWkm<-coordsW/1000
    #Transformation of coordinates in Long/Lat
    spW=SpatialPoints(coordsW)
    proj4string(spW)=proj4string(gridW)
    spWprojected=spTransform(spW, CRS("+init=epsg:4326"))
    lat_W=matrix(spWprojected@coords[,2])
    lon_W=matrix(spWprojected@coords[,1])
    coordsW=data.frame(long=lon_W,lat=lat_W)
    coordsW<-as.matrix(coordsW)
    acv <- autocov_dist(GD_mean, coordsW,type="inverse", nbs = neighbourhoodRadius,longlat=T)    
    return(acv)
}

## filter cells by number of individual, number of species
filter_by_ind_spec_taxcov <- function(Wr,thrshldInd,thrshldSpec) {
Wr.1=Wr[which(Wr$nb_indv >thrshldInd),]
### species
Wr.2=Wr.1[which(Wr.1$nb_species >(thrshldSpec-1)),]
    return(Wr.2)
}

## dataframe sensitivity results
len_df_sens=10
df_sens_empty=data.frame(filterSeq=c(0,1.9,1.9,1.9,1.9,3.9,1.9,8,4.9,2.9),
    filterSpec=c(0,2,3,8,2,2,3,2,2,2),
    numbCells=rep(NA,len_df_sens),
    propCells=rep(NA,len_df_sens),
    temperature=rep(NA,len_df_sens),
    slopeAvg=rep(NA,len_df_sens),
    regions=rep(NA,len_df_sens),
    autocor=rep(NA,len_df_sens),  
    ajustedRsquare=rep(NA,len_df_sens)
)

##########################################################################
## load data
dc = load_data("09-all_descripteurs/total_data_genetic_diversity_with_all_descripteurs.tsv")
dcT=dc[[1]];dcM=dc[[2]];dcF=dc[[3]];dcM.n1=dc[[4]];dcF.n1=dc[[5]];dcM.1=dc[[6]];dcF.1=dc[[7]]
## total number of species by cell
richness=read.csv("01-infos/equalarea_id_coordsCA_FWRS_MR_RS.csv")

##coeffs final model
## marine
modelMCoeffs=data.frame(temperature=0.341185,
regions=0.221148,
autocor=0.076247,
adj.r.squared=0.1598)
## freshwater
modelFCoeffs=data.frame(temperature=0.134435,
slopeAvg=-0.225852,
regions=0.121156,
autocor=0.237272,
adj.r.squared=0.1941)

##########################################################################
## load worldcoast shapefile
grid <- readOGR("01-infos/grid_equalarea200km","gridFish.b_260418")

##########################################################################
## MARINE
## add region indopacific/atlantic information for MARINE 
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
## format data
dcM.r2<-dcM.r1[,-c(5,8,9,10,11,13, 14,16,17,18,21,22,23)]
dcM.r3<-na.omit(dcM.r2)
## Adding distance to the coast as explanatory variable
coast<-read.csv("01-infos/distanceCote.csv", h=T,sep=",")
coast.r3=coast[which(coast$cell %in% dcM.r3$cell),]
dcM.r4<-data.frame(x=dcM.r3$x,
                  y=dcM.r3$y,cell=dcM.r3$cell,
                  GD_mean=dcM.r3$GD_mean,
                  nb_species=dcM.r3$nb_species,
                  nb_indv_mean=dcM.r3$nb_indv_mean,
                  cloMeanVal=dcM.r3$cloMeanVal,
                  bathyVal=dcM.r3$bathyVal,
                  O2=dcM.r3$O2,
                  SST=dcM.r3$SST,
                  regions=dcM.r3$regions,
                  distanceFromShore=coast.r3$Distance_from_shore)
### Missing value distance to the coast : 2 missing values
dcM.r4$distanceFromShore[c(139,337)]=mean(na.omit(dcM.r4)$distanceFromShore)
## add non normalised values
dcM.r4$nb_indv=dcM[which(dcM$cell %in% dcM.r4$cell),]$nb_indv_mean


##########################################################################
## FRESHWATER
## add region continent information for FRESHWATER
north_america=c("USA","MEX","CAN")
regions=countrycode(dcF.1$ISO3,"iso3c","continent")
regions[which(dcF.1$ISO3 %in% north_america)]="North_america"
regions[which(regions=="Americas")]="South_america"
regions[which(is.na(regions))]="Antarctica"
dcF.r1=cbind(dcF.1,regions)
names(dcF.r1)
#Put 0 for negative values of altitude
dcF.r1[which(dcF.r1$bathyVal < 0),]$bathyVal=0
## format data
dcF.r2= dcF.r1[,-c(5,8,9,10,11,12,13,14,16,17,18,19,21,22,23)] 
dcF.r3=na.omit(dcF.r2)   
dcF.r4=data.frame(GD_mean=c(dcF.r3$GD_mean),
                  nb_species=c(dcF.r3$nb_species),
                  nb_indv_mean=c(dcF.r3$nb_indv_mean),
                  SST=c(dcF.r3$SST),
                  Alt=c(scale(dcF.n1[which(dcF.n1$cell %in% dcF.r3$cell),]$bathyVal)),
                  regions=c(dcF.r3$regions),
                  cell=c(dcF.r3$cell),
                  x=c(dcF.r3$x),
                  y=c(dcF.r3$y))
## additional filters : Remove oceania and Antartica and bad anotated cells and run model again : 346 cells
dcF.r5=dcF.r4[dcF.r4$regions != 6,]
dcF.r5=dcF.r5[dcF.r5$regions != 2,]
dcF.r5=dcF.r5[dcF.r5$cell!=3304,]
dcF.r5=dcF.r5[dcF.r5$cell!=9146,]  #346 cells
## Other environmental variables (slopes, flow, ..) from  www.gebco.net
## Only slope range is significant. I keep only this variables
env<-read.csv("01-infos/EnvFreshwater.csv")
env1<-env[which(env$cell %in% dcF.r5$cell),]
env1$SlopeAvg[c(226,242,338)]=mean(na.omit(env1$SlopeAvg))
SlopeAvg=env1$SlopeAvg
dcF.r6<-cbind.data.frame(dcF.r5,SlopeAvg)
## add non normalised values
dcF.r6$nb_indv=dcF[which(dcF$cell %in%dcF.r6$cell),]$nb_indv_mean

##########################################################################
## MARINE: sensitivity analysis
sensitivityRes= df_sens_empty
sensitivityRes$neighbourhoodRadius=c(5000,5000,5000,7000,5000,7000,5000,7000,6000,5000)
for( i in seq(1,len_df_sens)) {
    # filter cells by #indv, #spec, taxonomic coverage    
    Mr.filterd=filter_by_ind_spec_taxcov(dcM.r4,
    sensitivityRes$filterSeq[i],
    sensitivityRes$filterSpec[i]
    )
    ## number of kept cells
    sensitivityRes$numbCells[i]=dim(Mr.filterd)[1]
    ## proportion of kept cells
    sensitivityRes$propCells[i]=sensitivityRes$numbCells[i]/dim(dcM.r4)[1]
    ## regression of autocorrelative variable against all the explicative variables
    acinv2a=autocorr_var(Mr.filterd,dcM,grid,sensitivityRes$neighbourhoodRadius[i] )
    modAutocovM=lm(acinv2a~SST+scale(bathyVal)+cloMeanVal+nb_species+nb_indv_mean+regions+scale(distanceFromShore),data= Mr.filterd)
    autocor=residuals(modAutocovM)
    modclimM.autof=lm(GD_mean~SST+regions+scale(autocor)+nb_species+nb_indv_mean,data= Mr.filterd)
    summaryMod=summary(modclimM.autof)
    ## SST
    sensitivityRes$temperature[i]=summaryMod$coefficients[2,1]
    sensitivityRes$temperaturePvalue[i]=summaryMod$coefficients[2,4]
    sensitivityRes$temperatureDiff[i]=distModCoeff(sensitivityRes$temperature[i], modelMCoeffs$temperature )
    ## regions
    sensitivityRes$regions[i]=summaryMod$coefficients[3,1]
    sensitivityRes$regionsPvalue[i]=summaryMod$coefficients[3,4]
    sensitivityRes$regionsDiff[i]=distModCoeff(sensitivityRes$regions[i], modelMCoeffs$regions)    
    ## autocor
    sensitivityRes$autocor[i]=summaryMod$coefficients[4,1]
    sensitivityRes$autocorPvalue[i]=summaryMod$coefficients[4,4]
    sensitivityRes$autocorDiff[i]=distModCoeff(sensitivityRes$autocor[i], modelMCoeffs$autocor)    
    ## adjusted R square (model performance)
    sensitivityRes$ajustedRsquare[i]=summaryMod$adj.r.squared
    sensitivityRes$ajustedRsquareDiff[i]=distModCoeff(summaryMod$adj.r.squared, modelMCoeffs$adj.r.squared)   
}
sensResM=t(sensitivityRes)

##########################################################################
## FRESHWATER : sensitivity analysis
sensitivityRes= df_sens_empty
sensitivityRes$neighbourhoodRadius=c(3000,3000,5000,5000,5000,5000,3000,5000,5000,3000)

for( i in seq(1,len_df_sens)) {
    # filter cells by #indv, #spec, taxonomic coverage    
    Fr.filterd=filter_by_ind_spec_taxcov(dcF.r6,
    sensitivityRes$filterSeq[i],
    sensitivityRes$filterSpec[i]
    )
    ## number of kept cells
    sensitivityRes$numbCells[i]=dim(Fr.filterd)[1]
    ## proportion of kept cells
    sensitivityRes$propCells[i]=sensitivityRes$numbCells[i]/dim(dcF.r6)[1]
    ## regression of autocorrelative variable against all the explicative variables
    acinv2a=autocorr_var(Fr.filterd,dcF,grid,sensitivityRes$neighbourhoodRadius[i] )
    modAutocovF=lm(acinv2a~SST+nb_species+nb_indv_mean+regions+Alt+SlopeAvg,data=Fr.filterd)
    autocor=residuals(modAutocovF)
    modclimF.autof=lm(GD_mean~SST+scale(SlopeAvg)+regions+scale(autocor)+nb_species+nb_indv_mean,data= Fr.filterd)
    summaryMod=summary(modclimF.autof)
    ## SST
    sensitivityRes$temperature[i]=summaryMod$coefficients[2,1]
    sensitivityRes$temperaturePvalue[i]=summaryMod$coefficients[2,4]
    sensitivityRes$temperatureDiff[i]=distModCoeff(sensitivityRes$temperature[i], modelFCoeffs$temperature)
    ## regions
    sensitivityRes$regions[i]=summaryMod$coefficients[4,1]
    sensitivityRes$regionsPvalue[i]=summaryMod$coefficients[4,4]
    sensitivityRes$regionsDiff[i]=distModCoeff(sensitivityRes$regions[i], modelFCoeffs$regions)
    ## autocor
    sensitivityRes$autocor[i]=summaryMod$coefficients[5,1]
    sensitivityRes$autocorPvalue[i]=summaryMod$coefficients[5,4]
    sensitivityRes$autocorDiff[i]=distModCoeff(sensitivityRes$autocor[i], modelFCoeffs$autocor)
    ## adjusted R square (model performance)
    sensitivityRes$ajustedRsquare[i]=summaryMod$adj.r.squared
    sensitivityRes$ajustedRsquareDiff[i]=distModCoeff(summaryMod$adj.r.squared, modelFCoeffs$adj.r.squared)
    ## slopeAvg
    sensitivityRes$slopeAvg[i]=summaryMod$coefficients[3,1]
    sensitivityRes$slopeAvgPvalue[i]=summaryMod$coefficients[3,4]
    sensitivityRes$slopeAvgDiff[i]=distModCoeff(sensitivityRes$slopeAvg[i], modelFCoeffs$slopeAvg)

}
sensResF=t(sensitivityRes)


##########################################################################
## write CSV files
sensResMround=round(sensResM,4)
write.csv(sensResMround,file="13-analysis/marine_sensitivity_model.csv",sep=",",col.names=F)
sensResFround=round(sensResF,4)
write.csv(sensResFround,file="13-analysis/freshwater_sensitivity_model.csv",sep=",",col.names=F)
