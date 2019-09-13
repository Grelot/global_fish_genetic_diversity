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
##
##
##########################################################################

library(sp)
library(spdep)

### from a data.frame (argument1)
### scale all the continuous explicative variables
### log function is applied on response variable "GD_mean"
### return a data.frame with formated data
format_dc <- function(dc) {
    dc.format= dc[which(dc$GD_mean != 0),]
    dc.format$HDI_2015=as.numeric(dc.format$HDI_2015)
    rownames(dc.format) <- NULL
    dc.format$AP=as.numeric(dc.format$AP)
    dc.format$x=scale(dc.format$x)
    dc.format$y=scale(dc.format$y)
    dc.format$VEL=scale(dc.format$VEL)
    dc.format$O2=scale(dc.format$O2)
    dc.format$SST=scale(dc.format$SST)
    dc.format$cloMeanVal=scale(dc.format$cloMeanVal)
    dc.format$nb_indv_mean=scale(dc.format$nb_indv_mean)
    dc.format$bathyVal=scale(dc.format$bathyVal)
    dc.format$fshD=scale(dc.format$fshD)
    dc.format$GD_mean=scale(log(dc.format$GD_mean))
    dc.format$GD_median=scale(log(dc.format$GD_median))
    dc.format$HDI_2015=scale(dc.format$HDI_2015)
    dc.format$cloMinVal=scale(dc.format$cloMinVal)
    dc.format$cloMaxVal=scale(dc.format$cloMaxVal)
    return(dc.format)
}

### from a linear model "mod"(argument1) and a name of file "nom_f"(argument2)
### it generates an I 's morgan figure into pdf format
I_morgan <- function(mod,nom_f) {
    xymod=xy[which(residuals(mod) != -9999),]
    voismod=dnearneigh(xymod,d1=0,d2=75,longlat=T,row.names=rownames(xymod))
    cormd=sp.correlogram(voismod,var=as.vector(residuals(mod)),order=4,method="I",randomisation=F,zero.policy=T)
    pdf(nom_f,width=15,height=10,paper='special')
    plot(cormd1)
    dev.off()
}

### load data from a .TSV file "tsvFile"(argument1)
### return 6 data.frames :
### cluster the data into 2 groups : marine and freshwater
### scale data (see function format_dc)
### remove cell with only one species
load_data <- function(tsvFile) {
    ##charger les donnees
    dcT<-read.table(tsvFile,header=T)
    #cluster marine
    dcM<-dcT[which(dcT$is_marine == 1),]
    #cluster freshwater
    dcF<-dcT[which(dcT$is_marine == 0),]
    #remove cells with only one species
    dcM.n1<-dcM[which(dcM$nb_species != 1),]
    dcF.n1<-dcF[which(dcF$nb_species != 1),]
    ##scale data
    dcM.1<-format_dc(dcM.n1)
    dcF.1<-format_dc(dcF.n1)
    return(list(dcT, dcM, dcF, dcM.n1, dcF.n1, dcM.1, dcF.1))
}


## build figure from a linear model of lines and points with ggplot
## it requires title and x lab and y lab arguments
## it requires linear model and x and y data
ggplot_lm<- function(mod,df_x,df_y,x_lab,y_lab,titre) {
    #R square
    rqs = summary(mod)$adj.r.squared
    #coefficients
    b = summary(mod)$coefficients[1]
    a1 = summary(mod)$coefficients[2]    
    df=data.frame(df_x,df_y)    
    dr_x=df_x
    dr_y=b+a1*dr_x
    dregr=data.frame(dr_x,dr_y)
    p= ggplot(df,aes(x=df_x,y=df_y))+
    theme_minimal()+
    geom_point(colour = "black", size = 2)+
    geom_line(data=dregr,aes(x=dr_x,y=dr_y), colour="red",size=1)+
    xlab(x_lab)+
    ylab(y_lab)+
    ggtitle(titre)
    #annotate("text",x=(max(df_x)/10)*8,y=(max(df_y)/10)*8, label = paste("y = ",round(b,digits=4),"+",round(a1,digits=4),"x\nr² = ",round(rqs,digits=4),"\n p-value < 0.001"))
    return(p)
}

## from shapefiles, it "plots" a world map
plot_world<- function(worldcoast,biglakes,riverlakes) {
  plot(rnorm(100),xlim=c(-167, 167),ylim=c(-84,75),xlab="",ylab="",type="n",axes=F)
  plot(worldcoast,add=TRUE,col="grey",border="black",lwd=0.2)
  plot(biglakes,add=TRUE,col="white",lwd=0.2)
  plot(riverlakes,add=TRUE,col="white",lwd=0.2)
  axis(side=1,line=0,cex.axis=1.2,lwd=0.35,tcl=-0.25,bg="white",labels=c("-150","-100","-50","0","50","100","150"), at=c(-150,-100,-50,0,50,100,150),mgp=c(3,0.5,0))
  axis(side=2,line=0,las=2,cex.axis=1.2,lwd=0.35,tcl=-0.25,bg="white",labels=c("-80","-60","-40","-20","0","20","40","60","80"), at=c(-80,-60,-40,-20,0,20,40,60,80),mgp=c(3, 0.5, 0))
}

## from a latitude coordinate, return 10-based latband value
coords_to_latband <- function(latcoord) {
    if(latcoord>=0) {
        bandcoord=10*floor(latcoord/10)+10
    } else{
        bandcoord=10*floor(latcoord/10)
    }
    return(bandcoord)
}

###############################################################################
## sensitiviy analysis functions

## distance between 2 model's coefficients
distModCoeff <- function(new,ref) {
    dist=(new-ref)/ref
    return(dist)
}

