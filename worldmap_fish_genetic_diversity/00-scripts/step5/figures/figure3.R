##########################################################################
## Codes for the paper:
## Analyse scripts for "Global patterns of fish genetic
## diversity increase with “current” temperature"
## Manel S, Guerin PE, Mouillot D, Velez L, Albouy C, Pellissier L
## 
## Montpellier, 2018, updated 2019
##
## functions used to generate figures
## Figure 3 : Determinant of the patterns of fish genetic diversity.
## Outputs of the linear model (lm) testing the effect of geographic,
## environmental, and sampling effect on the global pattern of
## marine (a-b) and freshwater (c-d) 
## Coefficients and confidence intervals for the factors of the models for 
## marine (a) and freshwater fishes (c). 
## Autocor is a spatial autocovariate that takes into account spatial
## autocorrelation in both our predicted and predictive variables. 
## Relative variance of genetic  diversity explained by
## the various factors estimated with the package hier.part in marine (b)
## and in freshwater fishes (d).
##
##########################################################################
## Libraries
lib_vect <-c("car","ggplot2","raster","rgdal","rgeos","MASS","hier.part","countrycode","lme4","sjPlot","plyr","gridExtra","visreg")
sapply(lib_vect,library,character.only=TRUE)


##########################################################################
## functions
source("00-scripts/step5/figures/functions.R")

convert_coordinates_meter_to_degree <- function(grid, coordsW) {
# convert meter-coordinates based on grid's projection into degree based on epsg:4326
    spW=SpatialPoints(coordsW)
    proj4string(spW)=proj4string(grid)
    spWprojected=spTransform(spW, CRS("+init=epsg:4326")) 
    lat_W=matrix(spWprojected@coords[,2])
    lon_W=matrix(spWprojected@coords[,1])
    degreeCoordsW=data.frame(long=lon_W,lat=lat_W) 
    return(degreeCoordsW)
}


##########################################################################
## load data
dc = load_data("09-all_descripteurs/total_data_genetic_diversity_with_all_descripteurs.tsv")
dcT=dc[[1]];dcM=dc[[2]];dcF=dc[[3]];dcM.n1=dc[[4]];dcF.n1=dc[[5]];dcM.1=dc[[6]];dcF.1=dc[[7]]


##########################################################################
## load worldcoast shapefile
grid <- readOGR("01-infos/grid_equalarea200km","gridFish.b_260418")


##########################################################################


##########################################################################
## FRESHWATER


##########################################################################
### add region continent information for FRESHWATER
north_america=c("USA","MEX","CAN")
regions=countrycode(dcF.1$ISO3,"iso3c","continent")
regions[which(dcF.1$ISO3 %in% north_america)]="North_america"
regions[which(regions=="Americas")]="South_america"
regions[which(is.na(regions))]="Antarctica"
dcF.r1=cbind(dcF.1,regions)
names(dcF.r1)
#Put 0 for negative values of altitude
dcF.r1[which(dcF.r1$bathyVal < 0),]$bathyVal=0
#cor.test(dcF.r1$SST,dcF.r1$bathyVal)


##########################################################################
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
write.table(dcF.r5,"09-all_descripteurs/models/dcF346Freshwater.txt",row.names=F)


##########################################################################
## Other environmental variables (slopes, flow, basin area..) from  www.gebco.net
## Only slope range is significant. I keep only this variables
env<-read.csv("01-infos/EnvFreshwater.csv")
env1<-env[which(env$cell %in% dcF.r5$cell),]
#mean(na.omit(env1$SlopeRange)) #641.68  NA: 226; 242 ; 338-->replace by mean values
env1$SlopeRange[c(226,242,338)]=642
SlopeRange=env1$SlopeRange
## missing data slope average mean=211 -->replace by mean values
env1$SlopeAvg[c(226,242,338)]=211
SlopeAvg=env1$SlopeAvg
env2<-cbind.data.frame(dcF.r5,cbind(SlopeRange,SlopeAvg))



##########################################################################
## Checking colinearity of variable with a VIF procedure in a simple LM
vif(lm(GD_mean~SST+nb_species+nb_indv_mean+regions+Alt+SlopeRange+SlopeAvg,data=env2))
## SlopeRange & SlopeAvg colinears
## The VIF is how much the variance of your regression coefficient is larger than it would otherwise have been if the variable had been completely uncorrelated with all the other variables in the model. 
vif(lm(GD_mean~SST+nb_species+nb_indv_mean+regions+Alt+SlopeAvg,data=env2))
#env2$SlopeRange=env2$SlopeAvg

##########################################################################
#Transformation of coordinates in Long/Lat
coordNS= data.frame(cbind(dcF$x,dcF$y,dcF$cell))
coordNS=coordNS[which(coordNS$X3 %in% dcF.r5$cell),]
coordsF<-data.frame(cbind(coordNS$X1, coordNS$X2)) #Coordinates in meter
coordsF=convert_coordinates_meter_to_degree(grid,coordsF)
coordsF <- as.matrix(coordsF)


##########################################################################
## regression of autocorrelative variable against all the explicative variables
acinv2a <- autocov_dist(dcF.r5$GD_mean,coordsF,type="inverse", nbs = 3000,longlat=T)
modAutocovF=lm(acinv2a~SST+nb_species+nb_indv_mean+regions+Alt+SlopeAvg,data= env2)
## residuals (part of spatial autocorrelation not explained by the model)
autocor=residuals(modAutocovF)


##########################################################################
## step AIC with term of spatial autocorrelation
modclimF.auto=lm(GD_mean~SST+nb_species+nb_indv_mean+regions+scale(autocor)+Alt+scale(SlopeAvg),data=env2)
stepAIC(modclimF.auto,trace = 1)  #default direction = backward
## final model (with spatial autocorrelation)
modclimF.autof=lm(GD_mean~SST+scale(SlopeAvg)+regions+scale(autocor)+nb_species+nb_indv_mean,data= env2)
summary(modclimF.autof)
extractAIC(modclimF.auto,k=2)
## check residuals
par(mfrow=c(2,2))
plot(modclimF.autof)
graphics.off()


##########################################################################
## check spatial autocorrelation of this new model
st.nb=knearneigh(coordsF,k=1,longlat = TRUE)
st.nb= knn2nb(knearneigh(coordsF, k = 1), row.names = rownames(coordsF))
st.lw=nb2listw(st.nb)
## Moran's test on residuals : non-parametric test of Moran Index (keep only this one at the end)
res.mod.ols <- moran.mc(as.matrix(residuals(modclimF.autof)), listw=st.lw, nsim=999)  
res.mod.ols # Of course non significant


##########################################################################
## confidence intervals
nom_variables=c("Mean regional\ntemperature","Average slope","Regions","Autocor","Number of species", "Number of sequences")
names(modclimF.autof$coefficients)=c("(Intercept)",nom_variables)
p_ic_F=plot_model(modclimF.autof,colors = c("blue","red"), type = "est",sort.est=T,vline.color = "grey", line.size=1.2,show.values = TRUE, value.offset = .25, value.size=8,digits=3)
p_ic_F=p_ic_F+theme(plot.title = element_text(size = 24),
                    axis.title=element_text(size=20),
                    axis.text=element_text(size=20),
                    plot.margin = unit(c(0.5,2,1,1), "cm"),
                    panel.grid.minor = element_blank(),
                    panel.grid.major = element_blank(),
                    panel.border=element_rect(fill=NA,size=1),
                    panel.background = element_blank())+ggtitle("(c)")


##########################################################################
## hier.part
nom_variables=c("Number of species", "Number of sequences","Mean regional\ntemperature","Regions","Average slope","Autocor")
dcF.r3.auto=cbind(env2,autocor)
hp_F=hier.part(dcF.r3.auto$GD_mean,dcF.r3.auto[,c(2,3,4,6,10,11)])
hp_df_F=data.frame(nom_var=nom_variables,hier_part=hp_F$I.perc$I)
hp_df_F$nom_var=factor(hp_df_F$nom_var,levels=hp_df_F$nom_var[order(hp_df_F$hier_part)])
p_hp_F=ggplot(data=hp_df_F, aes(x=nom_var, y=hier_part)) +
  geom_bar(stat="identity", fill="grey",colour="black")+
  ylim(0,100)+
  xlab("")+
  ylab("% Relative variance")+
  scale_y_discrete(limits=seq(0,100,20),expand = expand_scale(mult = c(0, .1)))+
  coord_cartesian(ylim = c(0, 100))+
  theme(plot.title = element_text(size = 24),
        axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=20,hjust=0.75,vjust=0.8,angle=45),
        axis.text.y = element_text(size=20),
        plot.margin = unit(c(0.5,1,1,1), "cm"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border=element_rect(fill=NA,size=1),
        panel.background = element_blank())+
  ggtitle("(d)")


##########################################################################
## MARINE


##########################################################################
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


##########################################################################
## format data
dcM.r2<-dcM.r1[,-c(5,8,9,10,11,13, 14,16,17,18,21,22,23)]
dcM.r3<-na.omit(dcM.r2)


##########################################################################
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


##########################################################################
## Selection top 200 cells
### Influence number of indviduals
length(which(dcM.r3$nb_indv_mean == 2)) #(>4.1 =more than 4 sequences for 146 cells (instead of 171))
dcM.r5<-dcM.r3[dcM.r3$nb_indv_mean > 2.9,]
write.table(dcM.r5, "09-all_descripteurs/models/S2I4C29M.txt",row.names=F)
### Influence of the number of species
length(which(dcM.r3$nb_species==2)) 
dcM.r5<-dcM.r3[dcM.r3$nb_species > 2.9,]
write.table(dcM.r5, "09-all_descripteurs/models/S8I2C36M.txt",row.names=F)


##########################################################################
## Checking colinearity of variable with a VIF procedure in a simple LM
vif(lm(GD_mean~SST+cloMeanVal+O2+scale(bathyVal)+nb_species+nb_indv_mean+regions+scale(distanceFromShore),data= dcM.r4)) 
vif(lm(GD_mean~SST+cloMeanVal+scale(bathyVal)+nb_species+nb_indv_mean+regions+scale(distanceFromShore),data= dcM.r4))


##########################################################################
#Transformation of coordinates in Long/Lat
coordNSM= data.frame(cbind(dcM$x,dcM$y,dcM$cell))
coordNSM=coordNSM[which(coordNSM$X3 %in% dcM.r4$cell),]
coordsM<-data.frame(cbind(coordNSM$X1, coordNSM$X2)) # Coordinates in meter
coordsMkm<-coordsM/1000 # Coordinates in kilometer
coordsM=convert_coordinates_meter_to_degree(grid,coordsM) # Coordinates in degree
coordsM <- as.matrix(coordsM)


##########################################################################
## regression of autocorrelative variable against all the explicative variables
acinv2a <- autocov_dist(dcM.r4$GD_mean,coordsM,type="inverse", nbs = 5000,longlat=T)
modAutocovM=lm(acinv2a~SST+bathyVal+cloMeanVal+nb_species+nb_indv_mean+regions+scale(distanceFromShore),data= dcM.r4) #O2 et SST ont VIF>30 : j'enleve O2; travailler sur .r2 ou .r3 revient au meme
## residuals (part of spatial autocorrelation not explained by the model)
autocor=residuals(modAutocovM)


##########################################################################
## step AIC with term of spatial autocorrelation
#modclimM.auto=lm(GD_mean~SST+nb_species+nb_indv_mean+regions+autocor,data=dcM.r4)
modclimM.auto=lm(GD_mean~SST+bathyVal+cloMeanVal+nb_species+nb_indv_mean+regions+scale(distanceFromShore)+scale(autocor),data=dcM.r4)
stepAIC(modclimM.auto)
### final model (with spatial autocorrelation)
modclimM.autof=lm(GD_mean~SST+regions+scale(autocor)+nb_species+nb_indv_mean,data= dcM.r4)
summary(modclimM.autof)
## check residuals
par(mfrow=c(2,2))
plot(modclimM.autof)
graphics.off()


##########################################################################
## check spatial autocorrelation of this model
st.nb.f=knearneigh(coordsM,k=1,longlat = TRUE)
st.nb.f= knn2nb(knearneigh(coordsM, k = 1), row.names = rownames(coordsM))
st.lw.f=nb2listw(st.nb.f)
res.mod.ols.f <- moran.mc(as.matrix(residuals(modclimM.autof)), listw=st.lw.f, nsim=999)
res.mod.ols.f


##########################################################################
## confidence intervals
nom_variables=c("Sea surface\ntemperature","Regions","Autocor","Number of species","Number of sequences")
names(modclimM.autof$coefficients)=c("(Intercept)",nom_variables)
p_ic_M=plot_model(modclimM.autof,colors = c("blue","red"), type = "est",sort.est=T,vline.color = "grey", line.size=1,show.values = TRUE, value.offset = .25, value.size=8,digits=3)
p_ic_M=p_ic_M+theme(plot.title = element_text(size = 24),
                    axis.title=element_text(size=20),
                    axis.text=element_text(size=20),
                    plot.margin = unit(c(0.5,2,1,1), "cm"),
                    panel.grid.minor = element_blank(),
                    panel.grid.major = element_blank(),
                    panel.border=element_rect(fill=NA,size=1),
                    panel.background = element_blank())+ggtitle("(a)")


##########################################################################
## hier.part
nom_variables=c("Sea surface\ntemperature","Regions","Autocor","Number of species","Number of sequences")
dcM.auto=cbind(dcM.r4,autocor)
hp_M=hier.part(dcM.auto$GD_mean,dcM.auto[,c(10,11,13,5,6)])
hp_df_M=data.frame(nom_var=nom_variables,hier_part=hp_M$I.perc$I)
hp_df_M$nom_var=factor(hp_df_M$nom_var,levels=hp_df_M$nom_var[order(hp_df_M$hier_part)])
p_hp_M=ggplot(data=hp_df_M, aes(x=nom_var, y=hier_part)) +
  geom_bar(stat="identity", fill="grey",colour="black")+
  ylim(0,100)+
  xlab("")+
  ylab("% Relative variance")+
  scale_y_discrete(limits=seq(0,100,20),expand = expand_scale(mult = c(0, .1)))+
  coord_cartesian(ylim = c(0, 100))+
  theme(plot.title = element_text(size = 24),
        axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=20,hjust=0.75,vjust=0.8,angle=45),
        axis.text.y = element_text(size=20),
        plot.margin = unit(c(0.5,1,1,1), "cm"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border=element_rect(fill=NA,size=1),
        panel.background = element_blank())+
  ggtitle("(b)")


##########################################################################


##########################################################################
## write pdf files
pdf("10-figures/figure3.pdf",width=16,height=14,paper='special')
grid.arrange(p_ic_M,p_hp_M,p_ic_F,p_hp_F, nrow=2, widths=c(3,4))
dev.off()
